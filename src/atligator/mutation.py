"""This module provides data structures and functions for prediction of mutations from pmatrices. These mutations can
be fed into external scoring programs after applying them to the scaffold structure.

:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-06
"""
import random
from copy import copy
from operator import itemgetter
from typing import List, Iterator, Tuple

from atligator.pdb_util import canonical_amino_acids
from atligator.prediction_matrix import BinderPredictionMatrix
from atligator.structure import CouplingInfo, DesignedBinder


class Mutation:
    """A mutation describes the exchange of a residue type at a defined location."""

    def __init__(self, residue_id: int, original_residue_type: str, mutated_residue_type: str, derived: bool = False,
                 representative: bool = False):
        """Creates a new instance.
        :param residue_id: the residue id within its chain
        :param original_residue_type: the original residue type to be changed by the mutation
        :param mutated_residue_type: the new residue type after mutation
        :param derived: whether this mutation has been derived from another mutation (e.g. by coupling)
        :param representative: whether this mutation is representative for a set of coupled mutations
        """
        self.residue_id = residue_id
        self.original_residue_type = original_residue_type
        self.mutated_residue_type = mutated_residue_type
        self.derived = derived
        self.representative = representative

    def can_combine(self, other: 'Mutation') -> bool:
        """This mutation can be combined with another mutation if both instances refer to different residue ids.
        :param other: the other mutation to check for combinability
        :return: whether the combination is valid
        """
        return not self.residue_id == other.residue_id

    def is_identity(self) -> bool:
        """
        :return: whether this mutation is identical, i.e., has no effect.
        """
        return self.original_residue_type == self.mutated_residue_type

    def get_coupled_mutations(self, coupling_info: List[CouplingInfo], designed_binder: DesignedBinder) \
            -> List['Mutation']:
        """
        Determines a list of mutations coupled with this mutation based on provided coupling information.
        :param coupling_info: describes which residues are coupled
        :param designed_binder: the region of the acomplex from which to take original residue types
        :return: a list of coupled mutations, including the context mutation itself
        """
        if coupling_info is not None:
            for ci in coupling_info:
                if len(ci[1]) < 2:
                    continue
                for resid in ci[1]:
                    if int(resid) == self.residue_id:
                        representant = ci[1][0]
                        original_restypes = ",".join(sorted(set([dr.original_residue_type for dr
                                                                 in designed_binder.residues
                                                                 if str(dr.residue_id) in ci[1]])))
                        return [Mutation(int(coupled_resid), original_restypes, self.mutated_residue_type,
                                         derived=(coupled_resid != representant),
                                         representative=(coupled_resid == representant))
                                for coupled_resid in ci[1]]
        return [self]

    def __str__(self) -> str:
        """converts a mutation into its string representation.
        :return: string serialization
        """
        return f"{self.original_residue_type}{self.residue_id}{'+' if self.representative else ''}" \
               f"{'/' if self.derived else ''}{self.mutated_residue_type}"

    def __eq__(self, other: 'Mutation') -> bool:
        """Considers this mutation's residue id, original type, and mutated type for comparison.
        :param other: the mutation to compare with
        :return: whether self and other are equal
        """
        return self.residue_id == other.residue_id and \
               self.mutated_residue_type == other.mutated_residue_type and \
               self.original_residue_type == other.original_residue_type and \
               self.derived == other.derived and self.representative == other.representative

    def equiv(self, other: 'Mutation') -> bool:
        """
        Performs an equality check ignoring the original residue type as well as the derived/representative flags.
        :param other: the mutation to compare with
        :return: wehther self and other are equivalent
        """
        return self.residue_id == other.residue_id and \
               self.mutated_residue_type == other.mutated_residue_type

    def __hash__(self) -> int:
        return hash(self.residue_id) + hash(self.original_residue_type) + hash(self.mutated_residue_type)


def read_mutation(word: str) -> Mutation:
    """Creates a mutation from its text serialization.
    :param word: the string to convert
    :return: a parsed runtime instance
    """
    residue_id = ""
    original_residue_type = ""
    mutated_residue_type = ""
    n_digits = 0
    derived = False
    coupled = False
    for c in word:
        if c.isdigit():
            n_digits += 1
            residue_id += c
        elif c == '+':
            coupled = True
        elif c == '/':
            derived = True
        elif n_digits < 1:
            original_residue_type += c
        else:
            mutated_residue_type += c
    return Mutation(int(residue_id), original_residue_type, mutated_residue_type,
                    representative=coupled, derived=derived)


class MutationSet:
    """A mutation set combines several mutations that do not contradict each other."""

    def __init__(self, mutations: List[Mutation]):
        """Creates a new mutation set.
        :param mutations: the initial content of the mutations list
        """
        self.mutations = mutations

    def combine(self, mutation: Mutation) -> 'MutationSet':
        """Combines this mutation set with a given mutation.
        :param mutation: the mutation to be combined with this mutation set
        :return: a combined mutation set, or None if the combination is illegal
        """
        if all(m.can_combine(mutation) for m in self.mutations):
            return MutationSet(self.mutations + [mutation])

    def crossbreed(self, other: 'MutationSet') -> 'MutationSet':
        """Combines the mutations of this mutation set with those of another one and returns a cross-breeded instance.
        In the case of conflicts, this mutation set is dominant, i.e., the mutated residue type of this mutation and
        the original residue type of the other mutation are accepted.
        :param other: the other mutation set to crossbreed with
        :return: a combined mutation set
        """
        muts = [copy(sm) for sm in self.mutations]
        if other is not None:
            for om in other.mutations:
                conflicting_muts = [sm for sm in self.mutations if sm.residue_id == om.residue_id]
                if len(conflicting_muts) == 0:
                    muts.append(copy(om))
                else:
                    for cm in conflicting_muts:
                        cm.original_residue_type = om.original_residue_type
        return MutationSet(muts)

    def __str__(self) -> str:
        """
        :return: a text serialization of this mutation set
        """
        result = ""
        for mut in sorted(self.mutations, key=lambda m: m.residue_id):
            if not mut.derived:
                result += f"{mut} "
        result = result[0:-1]
        return result

    def __eq__(self, other: 'MutationSet') -> bool:
        """Considers this mutation set's mutations for comparison.
        :param other: the mutation to compare with
        :return: whether self and other are equal
        """
        return len(self.mutations) == len(other.mutations) and all(mut in other.mutations for mut in self.mutations)

    def equiv(self, other: 'MutationSet') -> bool:
        """Considers this mutation set's mutations for equivalence. In contrast to equality comparison, the original
        residue types as well as the derived/representative flags of mutations are ignored.
        :param other: the mutation to compare with
        :return: whether self and other are equivalent
        """
        return len(self.mutations) == len(other.mutations) and (len(self.mutations) == 0 or all(
            any(m for m in other.mutations if m.equiv(mut)) for mut in self.mutations))

    def __hash__(self):
        return sum([hash(mut) for mut in self.mutations])

    def remove_mutation(self, res_id: int):
        for i_n, i in enumerate(self.mutations):
            if i.residue_id == res_id:
                del self.mutations[i_n]
                break


def read_mutation_set(line: str) -> MutationSet:
    """Recovers a mutation set from its string serialization.
    :param line: the string serialization
    :return: a mutation set containing the mutations parsed from this line
    """
    muts: List[Mutation] = []
    for word in line.split():
        muts.append(read_mutation(word))
    return MutationSet(muts)


def read_mutation_sets(filename: str) -> Iterator[MutationSet]:
    """Reads in a file that contains a collection of mutation sets, each represented by a line of text
    :param filename: the path of the file to open
    :return: a generator of mutation sets parsed from the file
    """
    with open(filename, 'r') as fo:
        for line in fo:
            yield read_mutation_set(line)


def recombine(mutation_sets: List[MutationSet], mutations: List[Mutation]) -> Iterator[MutationSet]:
    """Recombines a list of mutation sets with a list of base mutations. The resulting mutation set contains up to
    len(mutation_sets)*len(mutations) elements. All elements are valid mutations.
    :param mutation_sets: the mutation sets to recombine
    :param mutations: the mutations to apply to each mutation set if legal
    :return: a generator for recombined mutation sets. It does not yield the original contents of mutation_sets.
    """
    for ms in mutation_sets:
        for mut in mutations:
            combined_ms = ms.combine(mut)
            if combined_ms is not None:
                yield combined_ms


def find_best_mutations(pmatrix: BinderPredictionMatrix, n_mutations: int, allow_identities: bool = False) \
        -> List[Mutation]:
    """Finds the best n mutations suggested by a prediction matrix based on the quality of the mutation that can be
    deduced by each column. Only non-identical mutations are considered.
    :param pmatrix: the prediction matrix to analyze
    :param n_mutations: the number of mutations to find
    :param allow_identities: if true, pseudo-mutations that do not change resdue type are included in the result
    :return: a list of mutations (maximum size n_mutations) having best quality
    """
    mutations: List[Tuple[float, Mutation]] = []
    for line in pmatrix.lines:
        for col in line.columns:
            mut = Mutation(line.residue_id, line.original_restype, col.residue_type)
            if allow_identities or not mut.is_identity():
                mutations.append((col.quality(), mut))
    sorted_mutations = sorted(mutations, key=itemgetter(0), reverse=True)
    return list([e[1] for e in sorted_mutations])[0:n_mutations]


def generate_mutations(source: BinderPredictionMatrix or MutationSet, base_mutations: int, n_combinations: int,
                       allow_identities: bool = True) -> List[MutationSet]:
    """Generates a collection of mutation sets based on a binder prediction matrix.
    :param source: the prediction matrix to use as basis for mutations
    :param base_mutations: the number of atomic mutations to consider
    :param n_combinations: the maximum length of mutation sets created by recombination of base mutations
    :param allow_identities: if true, pseudo-mutations that do not change resdue type are included in the result
    :return: a list of mutation sets generated from this matrix and the specified parameters
    """
    if type(source) == BinderPredictionMatrix:
        mutations = find_best_mutations(source, base_mutations, allow_identities)
    else:
        mutations = list(mutation for mutation in source.mutations)[0:base_mutations]
    mutation_sets = list([MutationSet([m]) for m in mutations])
    for i in range(1, n_combinations):
        for ms in recombine(list(mutation_sets), mutations):
            if ms not in mutation_sets:
                mutation_sets.append(ms)
    return mutation_sets


# TODO Deprecated:
def select_random_matrix_mutation(pmatrix: BinderPredictionMatrix, allow_identities: bool = False) -> Mutation:
    """Selects a random mutation from this matrix. The probability for taking a specific random mutation is based on the
    quality factor of the respective column.
    :param pmatrix: taken as a basis for mutation selection
    :param allow_identities: if True, the selected mutation may be identical (e.g. SER123SER)
    :return: an instance of Mutation representing the random pick
    """

    def _random_pick() -> Mutation:
        """Utility function that performs a random pick in a pmatrix
        :return: the random generated mutation
        """
        quality_slice: float = random.uniform(0.0, pmatrix.total_quality())
        curr_qual: float = 0.0
        for pline in pmatrix.lines:
            for pcolumn in pline.columns:
                curr_qual += pcolumn.quality()
                if curr_qual > quality_slice:
                    return Mutation(int(pline.residue_id), pline.original_restype, pcolumn.residue_type)

    while True:
        mut = _random_pick()
        if mut is not None and (allow_identities or not mut.is_identity()):
            return mut


def select_random_mutation(pmatrix: BinderPredictionMatrix, allow_identities: bool = False,
                           allow_all: bool = False, random_mutations: bool = False, chance: float = 0.2) -> Mutation:
    """Selects a random mutation in either of two ways:
    By totally random selection of location and amino acid type OR by selection with pmatix probabilities.
    For totally random picks the probability for taking a specific random mutation is the same for
    every wild card and every mutation (except from PRO, CYS and GLY).
    The probability for taking a specific random mutation from the pmatrix is based on the
    quality factor of the respective column.
    :param pmatrix: Defines binder residues. Is also used for mutation selection.
    :param allow_identities: If True, the selected mutation may be identical (e.g. SER123SER)
    :param allow_all: Whether to allow mutations to all canonical amino acids. If False, Gly, Pro & Cys are not allowed.
    :param random_mutations: Whether to apply random mutations independent from pmatrix.
    :param chance: Chance of picking a totally random mutation
    :return: an instance of Mutation representing the random pick
    """

    def _random_pick() -> Mutation:
        """Utility function that performs a random pick
        :return: the random generated mutation
        """
        pline = random.choice(pmatrix.lines)
        return Mutation(int(pline.residue_id), pline.original_restype, random.choice(amino_acids))

    def _random_matrix_pick() -> Mutation:
        """Utility function that performs a random pick in a pmatrix
        :return: the random generated mutation
        """
        quality_slice: float = random.uniform(0.0, pmatrix.total_quality())
        curr_qual: float = 0.0
        for pline in pmatrix.lines:
            for pcolumn in pline.columns:
                curr_qual += pcolumn.quality()
                if curr_qual > quality_slice:
                    return Mutation(int(pline.residue_id), pline.original_restype, pcolumn.residue_type)

    if allow_all:
        amino_acids = canonical_amino_acids[:]
    else:
        amino_acids = canonical_amino_acids[:]
        amino_acids.remove("CYS")
        amino_acids.remove("PRO")
        amino_acids.remove("GLY")
    while True:
        if random_mutations and random.random() < chance:
            mut = _random_pick()
        else:
            mut = _random_matrix_pick()
        if mut is not None and (allow_identities or not mut.is_identity()):
            return mut
