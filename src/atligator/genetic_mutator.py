"""This module contains the implementation of a genetic mutation strategy that utilizes an atlas for calculating a
probability matrix for mutations dynamically, and a scorer for assessing the fitness of individuals.

The algorithm takes as input a list of suggested mutations. It applies these mutations to the scaffold structure in
isolation and determines whether the binding has improved based on the total energy calculated. Energy-increasing
mutations are then recombined for a given number of generations, where the individuals having the highest score
are given priority.

:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-28
"""

import logging as log
import os
import pickle
import random
from concurrent.futures import ThreadPoolExecutor
from copy import copy
from datetime import timedelta
from timeit import default_timer as timer
from typing import List, Tuple, Dict, Callable

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from atligator.acomplex import LigandBinderComplex
from atligator.acomplex import generate_complex, ComplexDescriptor
from atligator.atlas import Atlas
from atligator.base_scorer import BaseScorer
from atligator.deep_atlas import generate_deep_prediction_matrix, DeepAtlas
from atligator.mutation import MutationSet, select_random_mutation, Mutation
from atligator.pdb_util import mutate_residue, aa_3to1_conv
from atligator.prediction_matrix import generate_prediction_matrix, BinderPredictionMatrix, read_matrix
from atligator.similarity_scoring import generate_clash_pmatrix


class Individual:
    """Instances of this class represent invididuals of the population managed by a genetic mutator."""
    template_structure: str  # path to the PDB file this individual was derived from
    epoch: int  # the epoch in which the individual was created
    generation: int  # the generation within the epoch in which the individual was created
    iid: int  # a unique identifier of the individual that is usually managed by the mutator
    fitness: float  # The fitness calculated by the scorer, reflecting the fitness of the individual
    score_terms: Dict[str, float]  # score terms underlying the fitness
    optimized_structure: str  # PDB file that contains the optimized structure
    progenitor: 'Individual'  # the progenitor individual (belonging to the epoch)
    mother: 'Individual'  # mother (should exist for every non-scaffold individual)
    father: 'Individual'  # father (should exist for every non-scaffold and non-orphan individual)

    def __init__(self, mset: MutationSet):
        """Creates a new individual based on its mutation set
        :param mset: the set of mutations this individual reflects
        """
        self.mset = mset
        self.template_structure = None
        self.epoch = -1
        self.generation = -1
        self.iid = -1
        self.fitness = None
        self.score_terms = None
        self.optimized_structure = None
        self.progenitor = None
        self.mother = None
        self.father = None

    def crossbreed(self, father: 'Individual') -> 'Individual':
        """Combines two individuals by crossbreeding their mutation sets. Assumes the current individual as mother
        :param father: the individual to breed with
        :return: a cross-breeded individual
        """
        child = Individual(self.mset.crossbreed(father.mset))
        child.mother = self
        child.father = father
        return child

    def absolute_mset(self) -> MutationSet:
        """returns a mutation set that describes the mutations assigned to this individual, plus non-conflicting
        predecessor mutations from the progenitor individual, if any.
        :return: the absolute mutation set of this individual
        """
        if self.progenitor is None:
            return self.mset
        else:
            return self.mset.crossbreed(self.progenitor.absolute_mset())

    def __str__(self) -> str:
        return '{:10.4f} '.format(self.fitness) + '{:3} '.format(self.epoch) + \
               '{:3} '.format(self.generation) + '{:5} '.format(self.iid) + f"[{self.absolute_mset()}]" + " (" + \
               ' '.join([str(k) + ':' + '{:10.4f} '.format(v).strip() for k, v in self.score_terms.items()]) + ")"


class GeneticMutator:
    """A genetic mutator controls the generation of individuals over multiple epochs and generations. This class
    summarizes the data necessary for management and offers high-level controlling operations.
    """
    current_iid: int  # the id of the last individual generated
    generation: int  # keeps track of the current generation within the epoch
    population: List[Individual]  # the overall population
    msets: List[str]  # the overall population's msets in string representation for faster searching
    progens_by_epoch: List[Individual]  # list of progenitors; the index corresponds to the epoch in which used
    acomplex_by_epoch: List[LigandBinderComplex]  # acomplexes used for pmatrix derivation in every epoch
    pmatrix_by_epoch: List[BinderPredictionMatrix]  # pmatrices used for mutation prediction in every epoch
    previous_top_n: List[Individual]  # stores the best n_top individuals of the previous generation.

    def __init__(self, atlas: Atlas, deep_atlas: DeepAtlas, descriptor: ComplexDescriptor, n_top: int,
                 clustering_factor: float, clash_factor: float, deep_factor: float, random_mutations: bool,
                 random_mutation_chance: float, max_siblings: int, n_clones: int, delta_factor: float,
                 clash_merge_base: float, n_delta: float, std_dev_weight: float, n_quality_k: float):
        """Creates a new genetic mutator based on an atlas and user-defined input.
        :param atlas: the atlas to use for the creation of binder prediction matrices in every epoch
        :param deep_atlas: optional; adds deep learning facilities to pmatrixd calculation
        :param descriptor: user-defined input (ligand residue types, mutable host residues, couplings)
        :param n_top: the number of individuals to print after every epoch; a generation is to be terminated if this
        list does not change
        :param clustering_factor: probability of distance clustering based mutation selection
        :param clash_factor: probability to choose clash factoring for the selection of mutations
        :param deep_factor: probability to involve a deep learning strategy for mutation selection
        :param random_mutations: whether to apply random mutations.
        :param random_mutation_chance: The chance of random mutations (over ones based on pmatrix).
        :param max_siblings: the maximum number of individuals with identical mutation sets that may coexist per epoch
        :param n_clones: the number of clones automatically generated for new individuals (identical mutation sets).
        Original is also counted. -> n_clones 1 means individuals are not cloned!
        :param delta_factor: Correction factor for merging clash/quality with delta within calculation of mutation
        ranking: 'n_quality*base^-clash + n_delta/FACTOR'.
        :param clash_merge_base: Base of clash punishment for merging with normalized quality: n_quality*BASE^-clash.
        :param n_delta: Normalization constant for delta. Equals to delta value of 50%% in normalized state.
        :param std_dev_weight: Base in weighing the standard deviation: BASE^-std_dev * clash.
        :param n_quality_k: Normalization constant for quality. Equals to quality value of 50%% in normalized state.
        """
        self.atlas = atlas
        self.deep_atlas = deep_atlas
        self.descriptor = descriptor
        self.n_top = n_top
        self.clustering_factor = clustering_factor
        self.clash_factor = clash_factor
        self.deep_factor = deep_factor
        self.random_mutations = random_mutations
        self.random_mutation_chance = random_mutation_chance
        self.max_siblings = max_siblings
        self.n_clones = n_clones
        self.clash_params = clash_merge_base, delta_factor, n_delta, n_quality_k, std_dev_weight
        self.current_iid = 0
        self.generation = 0
        self.population = []
        self.msets = []
        self.progens_by_epoch = []
        self.acomplex_by_epoch = []
        self.pmatrix_by_epoch = []
        self.previous_top_n = []

    def current_epoch(self) -> int:
        """:return: the current epoch"""
        if len(self.progens_by_epoch) > 0:
            return len(self.progens_by_epoch)
        else:
            return 0

    def current_progenitor(self) -> Individual:
        """:return: the progenitor structure used in the current epoch"""
        if len(self.progens_by_epoch) > 0:
            return self.progens_by_epoch[-1]

    def current_acomplex(self) -> LigandBinderComplex:
        """:return: the acomplex used in the current epoch"""
        if len(self.acomplex_by_epoch) > 0:
            return self.acomplex_by_epoch[-1]

    def current_pmatrix(self) -> BinderPredictionMatrix:
        """:return: the pmatrix used in the current epoch"""
        if len(self.pmatrix_by_epoch) > 0:
            return self.pmatrix_by_epoch[-1]

    def initialize_deep_atlas(self):
        """Initializes the deep atlas of this mutator using the default parametration."""
        self.deep_atlas.train(self.atlas)

    def start_epoch(self, default_progenitor_path: str = None) -> None:
        """Starts a new epoch. Sets a new progenitor based upon the list of globally fittest individual such that the
        new progenitor does not have a mutation set indentical to the one of any preceding progenitor. It is assumed
        that update_acomplex and update_pmatrix are the first two operations to be called after starting a new epoch.
        :param default_progenitor_path: the structure to serve as progenitor in the new epoch in case there is no other
        progenitor available from the list of fittest individuals
        structure of all previous epochs is selected automatically
        :return the number of the current epoch
        """
        self.generation = 0
        fittest = self.sorted_population(False)
        progenitor = None
        for f in fittest:  # TODO hash?
            if not any(progen for progen in self.progens_by_epoch if progen.absolute_mset().equiv(f.absolute_mset())):
                progenitor = f
                break
        if progenitor is None:
            progenitor = Individual(MutationSet([]))
            progenitor.template_structure = default_progenitor_path
            progenitor.optimized_structure = default_progenitor_path
        progenitor.progenitor = self.current_progenitor()
        self.progens_by_epoch.append(progenitor)

    def update_acomplex(self, acomplexfile: str) -> LigandBinderComplex:
        """To be called explicitly after starting any epoch. Generates a new acomplex based on the current progenitor,
        saves it under the specified path, and returns it.
        :param acomplexfile: storage location of the newly created acomplex
        :return: the generated acomplex
        """
        progen_structure: Structure = PDBParser(QUIET=True).get_structure(self.current_progenitor().optimized_structure,
                                                                          self.current_progenitor().optimized_structure)
        acomplex = generate_complex(self.atlas, progen_structure, self.descriptor)
        with open(acomplexfile, 'wb') as fo:
            pickle.dump(acomplex, fo)
        self.acomplex_by_epoch.append(acomplex)
        return acomplex

    def update_pmatrix(self, pmatrixfile: str, **score_args) -> BinderPredictionMatrix:
        """To be called explicitly after updating the acomplex of any epoch. Generates a new pmatrix based on the
        current acomplex, saves it under the specified paths, and returns it.
        :param pmatrixfile: storage location of the newly generated acomplex
        :param score_args: will be passed to the scoring function
        :return: the generated pmatrix
        """
        pmatrix = generate_prediction_matrix(self.current_acomplex(), **score_args)
        pmatrix = pmatrix.normalize(criterion='quality', total=True) * self.clustering_factor
        if self.clash_factor > 0.0:
            clash_pmatrix = generate_clash_pmatrix(pmatrix, self.clash_params) \
                .normalize(criterion='quality', total=True)
            pmatrix = pmatrix.merge(clash_pmatrix * self.clash_factor)
        if self.deep_factor > 0.0:
            deep_pmatrix = generate_deep_prediction_matrix(self.current_acomplex(), self.deep_atlas) \
                .normalize(criterion='quality', total=True)
            pmatrix = pmatrix.merge(deep_pmatrix * self.deep_factor)
        pmatrix = pmatrix.normalize().sort_columns().sort_lines()
        if pmatrix.total_quality(True) > 0.0:
            pmatrix.save(pmatrixfile)
            self.pmatrix_by_epoch.append(pmatrix)
            return pmatrix

    def sorted_population(self, restrict_to_current_epoch: bool, max_elements: int = -1) -> List[Individual]:
        """Identifies the list of fittest individuals
        :param restrict_to_current_epoch: if True, individuals belonging to epochs different from the current are
        ignored
        :param max_elements: maximum length of the returned list
        :return: the fittest individuals with ascending score values
        """
        ce = self.current_epoch()
        spop = sorted([i for i in self.population if
                       (not restrict_to_current_epoch or i.epoch == ce) and (i.fitness is not None)],
                      key=lambda i: i.fitness)
        if max_elements > 0:
            return spop[:max_elements]
        else:
            return spop

    def start_generation(self) -> int:
        """Starts a new generation within the current epoch.
        :return: the generation number.
        """
        self.previous_top_n = self.sorted_population(True, self.n_top)
        self.generation += 1
        return self.generation

    def apply_random_mutations(self, ind: Individual, random_repeat: float) -> int:
        """Applies random mutations to the given individual until it fulfills sibling constraints and according to a
        given random repeat probability
        :param ind: the individual to apply random mutations to
        :param random_repeat: with this probability, random mutations are even applied if not necessary
        :return: the randomly mutated individual
        """
        n_rand = 0
        hashed_mset = mset_to_hash(ind.absolute_mset())
        while random.random() < random_repeat or \
                len([i for i in self.msets if i == hashed_mset]) >= self.max_siblings + 1 - self.n_clones:
            if len(ind.mset.mutations) > 0 and random.random() < 0.33:
                ran_mut: Mutation = random.choice(ind.mset.mutations)
                ran_cmuts = ran_mut.get_coupled_mutations(self.descriptor.coupling_info, self.current_acomplex().binder)
                for mut in ran_cmuts:
                    ind.mset.remove_mutation(mut.residue_id)
            else:
                rmut = select_random_mutation(self.current_pmatrix(), random_mutations=True,
                                              chance=self.random_mutation_chance)
                cmuts = rmut.get_coupled_mutations(self.descriptor.coupling_info, self.current_acomplex().binder)
                ind.mset = ind.mset.crossbreed(MutationSet(cmuts))
            n_rand += 1
            hashed_mset = mset_to_hash(ind.absolute_mset())
        return n_rand

    def create_progenitor_individual(self):
        """Generates a new individual that has exactly the properties of the progenitor.
        :return: the individual generated
        """
        ind = Individual(MutationSet([]))
        ind.progenitor = self.current_progenitor()
        ind.epoch = self.current_epoch()
        ind.generation = self.generation
        self.apply_random_mutations(ind, 0.0)
        return ind

    def create_orphan_individual(self, mother: Individual, random_repeat: float = 0.75) \
            -> Tuple[int, Individual]:
        """Generates a new individual based on an origin individual and a couple of random mutations.
        :param mother: the precursor individual
        :param random_repeat: probability of random mutations; repeatedly applied
        :return: a tuple consisting of the number of random mutations applied and the individual generated
        """
        if mother is None:
            mother = self.create_progenitor_individual()
        orphan = Individual(MutationSet(mother.mset.mutations[:]))
        orphan.mother = mother
        orphan.progenitor = self.current_progenitor()
        orphan.epoch = self.current_epoch()
        orphan.generation = self.generation
        n_rand = self.apply_random_mutations(orphan, random_repeat)
        return n_rand, orphan

    def create_breeded_individual(self, mother: Individual, father: Individual, random_repeat: float = 0.25) \
            -> Tuple[int, Individual]:
        """Generates a new individual by crossbreeding a given mother and father individual. If necessary or with a
        specified probability, random mutations are applied in addition. It is ensured that no identical individual
        is already managed by this mutator.
        :param mother: the mother individual; will be given priority in all mutation conflicts
        :param father: the father individual
        :param random_repeat: probability of random mutations; repeatedly applied
        :return: a tuple consisting of the number of random mutations applied and the individual breeded
        """
        child = mother.crossbreed(father)
        child.progenitor = self.current_progenitor()
        child.epoch = self.current_epoch()
        child.generation = self.generation
        n_rand = self.apply_random_mutations(child, random_repeat)
        return n_rand, child

    def has_improved_during_generation(self) -> bool:
        """:return: whether the top n individuals have not improved since finishing the last generation
        """
        return not self.previous_top_n == self.sorted_population(True)[:self.n_top]

    def insert_individual(self, ind: Individual) -> None:
        """Inserts an individual, assigns a new iid to it, and increments the current iid counter.
        Also stores a short version of total mutations set called "hash" in the mutator for faster checks.
        :param ind: the individual to insert
        """
        ind.iid = self.current_iid
        self.msets.append(mset_to_hash(ind.absolute_mset()))
        self.population.append(ind)
        self.current_iid += 1


def muts_to_hash(mutations: List[Tuple[int, str]]):
    hash = "" if len(mutations) >= 1 else "[]"
    mutation_list = []
    for mut in mutations:
        mutation_list.append((mut[0],
                              aa_3to1_conv[mut[1]]))
    mutation_list = sorted(mutation_list)
    for mut in mutation_list:
        hash += f"{mut[0]}{mut[1]}"
    return hash


def mset_to_hash(mset: MutationSet) -> str:
    hash = "" if len(mset.mutations) >= 1 else "[]"
    mutation_list = []
    for mut in mset.mutations:
        mutation_list.append((mut.residue_id,
                              aa_3to1_conv[mut.mutated_residue_type]))
    mutation_list = sorted(mutation_list)
    for mut in mutation_list:
        hash += f"{mut[0]}{mut[1]}"
    return hash


def mutate_mset(structure: Structure, chain_id: str, mset: MutationSet) -> None:
    """Applies one or more point mutations, described by a mutation set, to a given structure.
    :param structure: original structure not containing the mutation
    :param chain_id: the id of the chain where to apply the mutations
    :param mset: the mutations to apply
    """
    for mutation in mset.mutations:
        residue = structure[0][chain_id][int(mutation.residue_id)]
        mutate_residue(residue, mutation.mutated_residue_type)


def mutate_scaffold(scaffold_path: str, descriptor: ComplexDescriptor) -> None:
    """Applies the mutations described by the ligand info of a specified descriptor to a scaffold file. Saves the result
    under the same path.
    :param scaffold_path: a PDB file to which all ligand mutations are applied
    :param descriptor: its ligand info contains the ligand mutations to apply
    """
    parser = PDBParser(QUIET=True)
    scaffold_structure: Structure = parser.get_structure(scaffold_path, scaffold_path)
    for chain_id, residue_id, new_residue_type in descriptor.ligand_info:
        residue: Residue = scaffold_structure[0][chain_id][int(residue_id)]
        mutate_residue(residue, new_residue_type)
    io = PDBIO()
    io.set_structure(scaffold_structure)
    io.save(scaffold_path)


def genetic_mutate(scorer: BaseScorer, atlas: Atlas, deep_atlas: DeepAtlas or None, scaffold_path: str,
                   descriptor: ComplexDescriptor,
                   n_epochs: int, max_generations: int, generations_first_epoch: int, n_recomb: int, max_siblings: int,
                   n_top: int, n_clones: int, n_workers: int,
                   clustering_factor: float, score_threshold: float, distance_factor: float, orient_factor: float,
                   secor_factor: float, score_exponent: float,
                   clash_factor: float, clash_merge_base: float, delta_factor: float, n_delta: float,
                   std_dev_weight: float, n_quality_k: float,
                   deep_factor: float, random_mutations: bool, random_mutation_chance: float,
                   mutant_pattern: str, mutator_dumpfile: str, mutant_dir: str,
                   observer: Callable[[int, int, int], None] = None) -> GeneticMutator:
    """Performs an genetic mutation run based on a scaffold structure and a list of mutations to be applied to a
    dedicated host chain. If any of the mutations increases the total score, it will be taken into consideration for
    recombination with other successful mutations. This procedure is repeated until no more mutations are found or
    a given maximum number of generations has been reached. The output includes a ranking of all mutation sets attempted
    :param scorer: based on this instance of a subclass of BaseScorer the mutated structures are scored
    :param atlas: the Atlas to use for generating a new pmatrix in every epoch
    :param deep_atlas: the deep Atlas used to assist pmatrix calculation (optional)
    :param scaffold_path: the initial structure that does not contain any of the mutations
    :param descriptor: descriptor structure for ligand, binder, and coupled residues. From this descriptor, the
    protocol builds a new pmatrix and acomplex in every epoch
    :param n_epochs: the number of epochs (large iterations)
    :param max_generations: the maximum number of generations (small iterations) per epoch
    :param generations_first_epoch: the number of generations in first epoch. Equal max_generations if set 0
    :param n_recomb: the number of mutation sets to enqueue for recombination in each generation
    :param max_siblings: the maximum number of individuals with identical mutation sets that may coexist per epoch
    :param n_top: if the best n_top mutations do not change within a generation, the protocol terminates
    :param n_clones: this number of clones is created for every individual generated. Clones are then processed
    in parallel, each is carrying a unique individual ID. Original is also counted.
    :param n_workers: number of workers to use for parallel processing (-1 is for number of CPUs)
    :param clustering_factor: probability of distance clustering based mutation selection (default selection strategy)
    :param score_threshold: parameter for the scoring function: will be passed along while pmatrix calculation
    :param distance_factor: weight factor for distance coefficient of the scoring function. [1/Angstrom]
    :param orient_factor: weight factor for the primary orientation coefficient (C_alpha->C_beta) of the scoring
    function [1/radians]
    :param secor_factor: weight factor for the secondary orientation coefficient (C_alpha->C_O) of the scoring
    function [1/radians]
    :param score_exponent: the outcome of the scoring function is postprocessed by raising by this value
    :param clash_factor: probability to choose clash factoring for the selection of mutations in each selection round.
    :param clash_merge_base: Base of clash punishment for merging with normalized quality: n_quality*BASE^-clash.
    :param n_quality_k: Normalization constant for quality. Equals to quality value of 50%% in normalized state.
    :param n_delta: Normalization constant for delta. Equals to delta value of 50%% in normalized state.
    :param std_dev_weight: Base in weighing the standard deviation: BASE^-std_dev * clash.
    :param delta_factor: Correction factor for merging clash/quality with delta within calculation of mutation ranking:
    'n_quality*base^-clash + n_delta/FACTOR'.
    :param deep_factor: probability to involve a deep learning strategy for mutation selection
    :param random_mutations: whether to apply random mutations.
    :param random_mutation_chance: The chance of random mutations (over ones based on pmatrix).
    :param mutant_pattern: a file name pattern for mutated and scored structures. Must contain one $ symbol, which will
    be replaced by a unique identifier of the mutant
    :param mutator_dumpfile: if this file exists, the genetic mutator will be reconstructed from it. Furthermore, after
    every generation, the current progress is dumped into this file.
    :param mutant_dir: directory where mutant files are stored - relative to current directory.
    :param observer: callback function invoked after every recombination with current epoch, generation, recombination.
        None refers to return type.
    :return: a list of all individuals produced and scored in all epochs, ordered by score
    """
    def __dump_mutator(_mutator_to_dump: GeneticMutator) -> None:
        """Utility that pickles the mutator_to_dump.
        :param _mutator_to_dump: is pickled into mutator_dumpfile
        """
        if mutator_dumpfile is not None:
            with open(mutator_dumpfile, 'wb') as fo:
                pickle.dump(_mutator_to_dump, fo)

    def __load_mutator() -> GeneticMutator:
        """Utility that unpickles a mutator and returns it.
        Keep in mind that an epoch that has been discontinued will not be continued at the next generation
        after a restart. The next generation is epoch + 1 -> generation 1.
        :return: the instance unpickled from mutator_dumpfile
        """
        with open(mutator_dumpfile, 'rb') as fo:
            dumped_mutator = pickle.load(fo)
        return dumped_mutator

    def __mutate_and_score(_mut_indiv: Individual) -> None:
        """Nested worker function that applies the mutations described by one individual and scores the structure.
        :param _mut_indiv: the individual to mutate and score
        """
        job_start = timer()
        mutant_path = mutant_dir + ('' if mutant_dir.endswith('/') else '/') + \
                      mutant_pattern.replace("$", str(_mut_indiv.iid))
        log.info(f"\t\tCreating mutant {_mut_indiv.iid} based on template "
                 f"{_mut_indiv.template_structure}: ({_mut_indiv.mset})")
        job_parser = PDBParser(QUIET=True)
        template_structure: Structure = job_parser.get_structure(_mut_indiv.template_structure,
                                                                 _mut_indiv.template_structure)
        # TODO Make sure it works with more than one binder chain:
        mutate_mset(template_structure, descriptor.binder_info[0][0], _mut_indiv.absolute_mset())
        job_io = PDBIO()
        job_io.set_structure(template_structure)
        job_io.save(mutant_path)
        log.info(f"\t\tScoring mutant {_mut_indiv.iid}...")
        mutant_score = scorer.score(mutant_path, descriptor)
        if os.path.isfile(mutant_path):  # TODO make sure unrelaxed file is not needed any more
            os.remove(mutant_path)
        else:
            log.warning(f"Error: {mutant_path} file not found")
        job_end = timer()
        log.info(f"\t\tDone with mutant {_mut_indiv.iid}. Time elapsed: "
                 f"{str(timedelta(seconds=job_end - job_start)).split('.')[0]}. Total score: {mutant_score}")
        _mut_indiv.optimized_structure = mutant_score.optimized_structure
        _mut_indiv.fitness = mutant_score.score
        _mut_indiv.score_terms = mutant_score.properties

    if not os.path.exists(mutant_dir):
        os.mkdir(mutant_dir)
    if n_workers < 0:
        n_workers = os.cpu_count()
    if generations_first_epoch < 0:
        generations_first_epoch = max_generations
    log.info("Initialization of genetic protocol.")
    mutator: GeneticMutator = None
    try:
        mutator = __load_mutator()
        log.info(f"\tRecovered dumped mutator instance: {mutator_dumpfile}. Starting resume mode overriding parameters "
                 f"as requested.")
        mutator.atlas = atlas
        mutator.deep_atlas = deep_atlas
        mutator.descriptor = descriptor
        mutator.n_top = n_top
        mutator.clash_factor = clash_factor
        mutator.deep_factor = deep_factor
        mutator.random_mutations = random_mutations
        mutator.random_mutation_chance = random_mutation_chance
        mutator.max_siblings = max_siblings
        mutator.n_clones = n_clones
        mutator.clash_params = clash_merge_base, delta_factor, n_delta, n_quality_k, std_dev_weight
    except FileNotFoundError:
        log.info(f"\tNo such dumpfile found: {mutator_dumpfile}. Skipping resume mode.")
    except EOFError:
        log.error(f"\tError while unpickling dumped mutator {mutator_dumpfile}. Disabling dumping for this run!")
        mutator_dumpfile = None
    except TypeError:
        pass
    if mutator is None:
        log.info(f"\tCreating a new mutator instance based on parameters passed.")
        mutator = GeneticMutator(atlas, deep_atlas, descriptor, n_top, clustering_factor, clash_factor, deep_factor,
                                 random_mutations, random_mutation_chance, max_siblings, n_clones, clash_merge_base,
                                 delta_factor, n_delta,
                                 std_dev_weight, n_quality_k)
        log.info(f"\tApplying ligand mutations to scaffold structure {scaffold_path}...")
        mutate_scaffold(scaffold_path, descriptor)
    if deep_factor > 0.0:
        log.info("\tInitializing deep atlas...")
        mutator.initialize_deep_atlas()
    log.info("Initialization done.")
    if mutator.current_epoch() > 0 or mutator.generation > 0:
        log.info(f"Resuming dumped protocol {mutator_dumpfile} after {mutator.current_epoch()} epochs, "
                 f"{mutator.generation} generations.")
        unscored_indivs = [i for i in mutator.population if i.fitness is None]
        if len(unscored_indivs) > 0:
            log.info(f"There are {len(unscored_indivs)} individuals left unscored. Catching up on this first.")
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                for _ in executor.map(__mutate_and_score, unscored_indivs):
                    pass
    while mutator.current_epoch() < n_epochs:
        mutator.start_epoch(scaffold_path)
        log.info(f"Starting epoch {mutator.current_epoch()} of {n_epochs} "
                 f"using progenitor structure {mutator.current_progenitor().optimized_structure}.")
        log.info("Creating acomplex...")
        acomplexfile: str = f"epoch{mutator.current_epoch()}.acomplex"
        mutator.update_acomplex(acomplexfile)
        log.info(f"Stored acomplex under {acomplexfile}.")
        log.info("Creating pmatrix...")
        pmatrixfile: str = f"epoch{mutator.current_epoch()}.pmatrix"
        if os.path.isfile(pmatrixfile):
            mutator.pmatrix_by_epoch.append(read_matrix(pmatrixfile))
        else:
            if mutator.update_pmatrix(pmatrixfile,
                                      score_threshold=score_threshold,
                                      distance_factor=distance_factor,
                                      orient_factor=orient_factor,
                                      secor_factor=secor_factor,
                                      score_exponent=score_exponent) is None:
                log.warning(f"Generated pmatrix shows a total quality of 0.0. Aborting epoch {mutator.current_epoch()} "
                            f"after 0 generations.")
                break
        log.info(f"Stored pmatrix under {pmatrixfile}.")
        while mutator.generation < (generations_first_epoch if mutator.current_epoch() == 1 else max_generations):
            mutator.start_generation()
            log.info(f"\tGeneration {mutator.generation} of "
                     f"{generations_first_epoch if mutator.current_epoch() == 1 else max_generations} in epoch "
                     f"{mutator.current_epoch()} of {n_epochs}.")
            recomb_indivs: List[Individual] = []
            for recomb in range(0, n_recomb):
                log.info(f"\t\tRecombination {recomb + 1} of {n_recomb}.")
                fittest = mutator.sorted_population(True)
                template = mutator.current_progenitor()
                progen_indiv = mutator.create_progenitor_individual()
                n_rand = 0
                if mutator.current_iid == 0:
                    log.info("\t\t\tFilling up with scaffold.")
                    indiv = progen_indiv
                elif recomb >= len(fittest) or recomb >= min(n_recomb - 1, int(n_recomb * 0.75)):
                    log.info("\t\t\tFilling up with orphan child.")
                    n_rand, indiv = mutator.create_orphan_individual(progen_indiv, 0.75)
                else:
                    mother = fittest[recomb]
                    if random.random() < 0.25:
                        father = progen_indiv
                    else:
                        father_index = int(random.uniform(0, min(len(fittest),
                                                                 max(int(len(fittest) / 4), 4 * n_recomb))))
                        father = fittest[father_index]
                    log.info(f"\t\t\tMother (score={mother.fitness}): [{mother.mset}]")
                    log.info(f"\t\t\tFather (score={father.fitness}): [{father.mset}]")
                    n_rand, indiv = mutator.create_breeded_individual(mother, father, 0.25)
                    template = sorted([progen_indiv, mother, father],
                                      key=lambda i: float('inf') if i.fitness is None else i.fitness)[0]
                log.info(f"\t\t\tChild after {n_rand} random mutation(s): [{indiv.mset}]")
                indiv.template_structure = template.optimized_structure if template.optimized_structure is not None \
                    else template.template_structure
                mutator.insert_individual(indiv)
                recomb_indivs.append(indiv)
                for i in range(1, n_clones):
                    clone = copy(indiv)
                    mutator.insert_individual(clone)
                    log.info(f"\t\tCreated individual {clone.iid} as a clone of {indiv.iid}.")
                    recomb_indivs.append(clone)
            n = len(recomb_indivs)
            log.info(f"\tParallel processing of {n} mutants...")
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                i = 0
                for _ in executor.map(__mutate_and_score, recomb_indivs):
                    i += 1
                    if observer is not None:
                        observer(mutator.current_epoch(), mutator.generation, i)
            log.info("\tParallel processing done.")
            __dump_mutator(mutator)
            top_n_pop = mutator.sorted_population(True, mutator.n_top)
            log.info(f"\tTop {len(top_n_pop)} population after epoch {mutator.current_epoch()}, "
                     f"generation {mutator.generation}:")
            progen_mset = mutator.current_progenitor().absolute_mset()
            if len(progen_mset.mutations) > 0:
                log.info(f"\tAll individuals inherit the progenitor mutation set: [{progen_mset}]")
            for indiv in top_n_pop:
                log.info('\t\t{:10.4f} '.format(indiv.fitness) + '{:3} '.format(indiv.generation) +
                         '{:5} '.format(indiv.iid) + f"[{indiv.mset}]")
            if mutator.generation < (generations_first_epoch if mutator.current_epoch() == 1 else max_generations) \
                    and not mutator.has_improved_during_generation():
                log.warning(f"No improvement over previous top {n_top} ranking. Aborting epoch "
                            f"{mutator.current_epoch()} after {mutator.generation} generations.")
                break
        ov_top_n_pop = mutator.sorted_population(False, mutator.n_top)
        log.info(f"Overall top {len(ov_top_n_pop)} population after epoch {mutator.current_epoch()}:")
        for indiv in ov_top_n_pop:
            log.info(f"\t{indiv}")
    log.info("Genetic protocol done.")
    return mutator


class NegativeIndividual(Individual):
    """Instances of this class represent invididuals of the population for negative design."""

    def __init__(self, mset: MutationSet, nd_mset: MutationSet):
        """Creates a new individual based on its mutation set
        :param mset: the set of binder mutations this individual reflects
        :param nd_mset: the set of ligand mutations this individual reflects for negative design approach
        """
        super().__init__(mset)
        self.nd_mset = nd_mset
        self.nd_iid = None
        # other score sets are scores that result in repetition of scoring. Including fitness and score terms.
        self.other_score_sets: List[Tuple[float, Dict[str, float]]] = []

    def __str__(self) -> str:
        # return '{:10.4f} '.format(self.fitness) + '{:5} '.format(self.iid) + f"[{self.nd_mset}]" + " (" + \
        #        ' '.join([str(k) + ':' + '{:10.4f} '.format(v).strip() for k, v in self.score_terms.items()]) + ")"
        first = f"    [{self.nd_mset}] " + '{:10.4f} '.format(self.fitness) + " (" + \
               ' '.join([str(k) + ':' + '{:10.4f} '.format(v).strip() for k, v in self.score_terms.items()]) + ")"
        other = ""
        for score_set in self.other_score_sets:
            other += "\n" + " "*10 + '{:10.4f} '.format(score_set[0]) + " (" + \
                     ' '.join([str(k) + ':' + '{:10.4f} '.format(v).strip() for k, v in score_set[1].items()]) + ")"

        return first + other
