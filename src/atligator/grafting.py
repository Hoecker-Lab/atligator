"""This module contains classes and functions for the automatic or manual grafting of pockets based on data obtained by
the pocket miner and user-defined input (ligand binder complex descriptors). When several residue pockets are grafted,
the algorithms take care that the globally optimal graft is preferred. All in all, the methods implemented here are very
coarse-grained and the PDB structures generated are physically improper and demand for manual revision and/or
computational optimization.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-06-11
"""
import copy
from itertools import product, combinations
from operator import itemgetter
from typing import Callable, Union
from typing import List, Tuple, Dict, Set

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from atligator.acomplex import LigandBinderComplex, ComplexDescriptor, generate_complex
from atligator.atlas import AtlasDatapoint
from atligator.pocket_miner import Pocket, PocketCluster
from atligator.structure import DesignedBinderResidue, get_icoor_from_aresidue, DesignedLigandResidue, AbstractResidue


class ResidueGraft:
    """
    A residue graft describes the change of a residue type and the coordinates of its atoms, based on datapoints from
    a pocket cluster.
    """

    def __init__(self, binder_res: DesignedBinderResidue, pocket_cluster: PocketCluster, penalty: float):
        """Creates a new instance.
        :param binder_res: the binder residue to be mutated from the acomplex
        :param pocket_cluster: the pocket cluster that provides type and atom coordinates of the grafted residue
        :param penalty: deviation of the cluster centroid from the binder residue's position and orientation
        """
        self.binder_res = binder_res
        self.pocket_cluster = pocket_cluster
        self.penalty = penalty

    def __repr__(self):
        return f"<ResidueGraft {self.binder_res.original_residue_type}{self.binder_res.residue_id}: " \
               f"{self.pocket_cluster.binder_restype}>"


class PocketGraft:
    """
    A pocket graft represents a set of binder residue grafts, which have been determined based on the properties of
    one given ligand residue.
    """

    def __init__(self, ligand_res: DesignedLigandResidue, residue_grafts: List[ResidueGraft]):
        """Creates a new instance.
        :param ligand_res: the ligand residue whose surrounding binder residues are subject to grafting
        :param residue_grafts: the residue grafts contained by this pocket graft.
        """
        self.ligand_res = ligand_res
        self.residue_grafts = residue_grafts

    def __repr__(self):
        return f"<PocketGraft {self.ligand_res.restype}{self.ligand_res.residue_id}: {self.residue_grafts}>"

    def penalty(self) -> float:
        """
        :return: Mean penalty of all residue grafts contained by this pocket graft
        """
        if len(self.residue_grafts) == 0:
            return 0.0
        return sum([rg.penalty for rg in self.residue_grafts]) / float(len(self.residue_grafts))


class BinderGraft:
    """
    A binder graft consists of several residue grafts, grouped by pocket grafts.
    """

    def __init__(self, pocket_grafts: List[PocketGraft]):
        """Creates a new instance
        :param pocket_grafts: disjoint partition of residue grafts.
        """
        self.pocket_grafts = pocket_grafts

    def __repr__(self):
        return f"<BinderGraft {self.pocket_grafts}>"

    def penalty(self) -> float:
        """
        :return: Mean penalty of all pocket grafts contained by this binder graft
        """
        if len(self.pocket_grafts) == 0:
            return 0.0
        return sum([pg.penalty() for pg in self.pocket_grafts]) / float(len(self.pocket_grafts))


def calculate_penalty(binder_res: DesignedBinderResidue, cluster_centroid: AtlasDatapoint,
                      distance_factor: float, orient_factor: float, secor_factor: float) -> float:
    """Calculates the distance between a wildcard residue and the centroid of a cluster.
    :param binder_res:
    :param cluster_centroid:
    :param distance_factor: weight factor for the distance term [1/Angsrom]
    :param orient_factor: weight factor for primary orientation [1/radians]
    :param secor_factor: weight factor for secondary orientation [1/radians]
    :return: a positive real value reflecting the penalty
    """
    calpha_distance: float = abs((binder_res.calpha() - cluster_centroid.binder_calpha()).norm())
    orient_angle: float = abs(binder_res.orientation().angle(cluster_centroid.binder_orientation()))
    sec_orient_angle: float = abs(binder_res.secondary_orientation().angle(
        cluster_centroid.binder_secondary_orientation()))
    return distance_factor * calpha_distance + orient_factor * orient_angle + secor_factor * sec_orient_angle


def get_optimal_pocket_graft(ligand_residue: DesignedLigandResidue, wildcard_residues: List[DesignedBinderResidue],
                             pocket: Pocket, penalty_threshold: float,
                             distance_factor: float, orient_factor: float, secor_factor: float) -> PocketGraft:
    """Yields the optimal assignment between wildcard residues and pocket members, such that the sum of penalties
    between the wildcard residue and the pocket's centroid is globally minimal. The returned assignment in terms of
    pocket grafts may leave wildcard residues unbound
    :param ligand_residue: the ligand residue where the pocket is to be grafted
    :param wildcard_residues: to be assigned to pocket members
    :param pocket: its members are assigned wildcard residues
    :param penalty_threshold: if the penalty exceeds this threshold, corresponding pocket grafts are ignored
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    :return: a pocket graft that reflects the calculated assignment
    """
    penalties: Dict[Tuple[int, int], float] = {}
    for i, wcres in enumerate(wildcard_residues):
        for j, cluster in enumerate(pocket.clusters):
            pen = calculate_penalty(wcres, cluster.centroid, distance_factor, orient_factor, secor_factor)
            penalties[(i, j)] = pen
    occupied_is: Set[int] = set()
    occupied_js: Set[int] = set()
    residue_grafts: List[ResidueGraft] = []
    for (i, j), pen in sorted(list(penalties.items()), key=itemgetter(1)):
        if len(residue_grafts) == len(wildcard_residues) or pen > penalty_threshold:
            break
        if i not in occupied_is and j not in occupied_js:
            residue_graft = ResidueGraft(wildcard_residues[i], pocket.clusters[j], pen)
            residue_grafts.append(residue_graft)
            occupied_is.add(i)
            occupied_js.add(j)
    return PocketGraft(ligand_residue, residue_grafts)


def get_pocket_grafts_for_candidates(ligand_residue: DesignedLigandResidue,
                                     wildcard_residues: List[DesignedBinderResidue], pocket_candidates: List[Pocket],
                                     penalty_threshold: float, distance_factor: float, orient_factor: float,
                                     secor_factor: float) -> List[PocketGraft]:
    """Calculates the optimal pocket graft for every element in a list of candidate pockets. Orders the candidates by
    the overall penalty their grafting would imply (assuming that nothing has been grafted yet)
    :param ligand_residue: see get_optimal_pocket_graft
    :param wildcard_residues: see get_optimal_pocket_graft
    :param pocket_candidates: candidates whose optimal grafts are calculated and ordered by penalty
    :param penalty_threshold: see get_optimal_pocket_graft
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    :return: a list of pocket grafts corresponding to the candidates passed, ordered by increasing penalty
    """
    return sorted([get_optimal_pocket_graft(ligand_residue, wildcard_residues, p, penalty_threshold, distance_factor,
                                            orient_factor, secor_factor) for p in pocket_candidates],
                  key=lambda pg: pg.penalty())


def get_optimal_binder_graft(pocket_grafts: Dict[int, List[PocketGraft]], max_binder_grafts: int = 1) \
        -> List[BinderGraft]:
    """ Given a dictionary of containing lists of pocket grafts, which potentially contain contradicting grafts
    for equivalent residues, this function selects the binder grafts with minimal penalty containing as many as possible
    pocket grafts (one per ligand position).
    :param pocket_grafts: pocket grafts obtained by local optimization of ligand pockets. In a dictionary with indices
    (range(x)) as keys (for the different ligand positions) and a list of pocket grafts as values (ordered by penalty
    - should not be required though).
    :param max_binder_grafts: Maximum number of binder grafts returned
    :return: a list of minimum-penalty and non-contradictory binder grafts, sorted by ascending penalty
    """

    def get_combinations(nones: int = 0):
        if not nones:
            return product(*(range(len(x)) for x in pocket_grafts.values()))

        none_combinations = combinations(range(n_pockets), nones)
        combin = []
        for none_positions in none_combinations:
            x_list = list(range(len(x)) for k, x in pocket_grafts.items() if k not in none_positions)
            for none in none_positions:
                x_list.insert(none, [None])
            combin.extend(list(product(*x_list)))
        return combin

    def get_penalty(graft_ids):
        return sum(penalty_dict[gi_n][gi] for gi_n, gi in enumerate(graft_ids))

    def mutually_exclusive(c):
        sets = [residue_sets[lg][pg] for lg, pg in enumerate(c)]
        return not any(set_c[0] & set_c[1] for set_c in combinations(sets, 2))

    residue_sets = {}
    penalty_dict = {}
    for lig, pgs in pocket_grafts.items():
        residue_sets[lig] = {}
        penalty_dict[lig] = {}
        for i, pg in enumerate(pgs):
            residue_sets[lig][i] = {(rg.binder_res.chain, rg.binder_res.residue_id) for rg in pg.residue_grafts}
            penalty_dict[lig][i] = pg.penalty()
        residue_sets[lig][None] = set()
        penalty_dict[lig][None] = 0.0

    n_pockets = len(list(pocket_grafts))
    binder_grafts = []
    n_lacking_pockets = 0
    while n_lacking_pockets < n_pockets:
        for comb in sorted(get_combinations(n_lacking_pockets), key=get_penalty):
            if mutually_exclusive(comb):
                binder_grafts.append(
                    BinderGraft(pocket_grafts=[pocket_grafts[ci][c] for ci, c in enumerate(comb) if c is not None])
                )
                if len(binder_grafts) >= max_binder_grafts:
                    return binder_grafts
        n_lacking_pockets += 1
    return binder_grafts


def select_best_pocket_graft(candidates: List[PocketGraft]) -> PocketGraft:
    """Default strategy for selecting an element from the list of pocket candidates (ordered by penalty). Take the
    cheapest element.
    :param candidates:
    :return: the first (i.e., lowest-penalty) candidate
    """
    return candidates[0]


def get_optimal_graft_for_acomplex(acomplex: LigandBinderComplex, pockets: Dict[str, List[Pocket]],
                                   penalty_threshold: float, distance_factor: float, orient_factor: float,
                                   secor_factor: float,
                                   select: Callable[[List[BinderGraft]], BinderGraft] = select_best_pocket_graft,
                                   loc_observer: Callable[[int, int], None] = None,
                                   glob_observer: Callable[[int, int], None] = None) -> BinderGraft or None:
    """ Based on a ligand/binder complex and a pocket dictionary, this method calculates an optimal binder graft based
    on a configurable candidate selection strategy
    :param acomplex: the ligand binder complex that defines both ligand residues and wildcards in the binder
    :param pockets: the output of a pocket mining run providing the basis for pocket selection
    :param penalty_threshold: see get_optimal_pocket_graft
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    :param select: the selection strategy to apply for choosing a binder graft
    :param loc_observer: callback to react to the progress of the local optimization stage
    :param glob_observer: callback to react to the progress of the local optimization stage
    :return: a globally optimal binder graft, or None if selection was cancelled
    """
    pocket_grafts: Dict[int, List[PocketGraft]] = {}
    loc_len = len(acomplex.ligand.residues)
    if loc_observer is not None:
        loc_observer(0, loc_len)
    for i, ligres in enumerate(acomplex.ligand.residues):
        icoor = get_icoor_from_aresidue(ligres)
        wc_res = [br.transform_coor(icoor, True) for br in acomplex.binder.residues]
        ligpockets = pockets[ligres.restype]
        pgs = get_pocket_grafts_for_candidates(ligres, wc_res, ligpockets, penalty_threshold, distance_factor,
                                               orient_factor, secor_factor)
        pocket_grafts[i] = pgs
        if loc_observer is not None:
            loc_observer(i, loc_len)

    binder_grafts = get_optimal_binder_graft(pocket_grafts=pocket_grafts, max_binder_grafts=1)
    binder_graft = select(binder_grafts)
    if binder_graft is None:
        return None
    glob_len = len(binder_graft.pocket_grafts)
    if glob_observer is not None:
        glob_observer(0, glob_len)
    for i, pg in enumerate(binder_graft.pocket_grafts):
        icoor = get_icoor_from_aresidue(pg.ligand_res)
        for rg in pg.residue_grafts:
            binder_res: DesignedBinderResidue = rg.binder_res
            rg.binder_res = binder_res.transform_coor(icoor, False)
        if glob_observer is not None:
            glob_observer(i, glob_len)
    return binder_graft


def apply_residue_graft(binder_residue: Residue, residue_graft: ResidueGraft,
                        distance_factor: float, orient_factor: float, secor_factor: float) -> None:
    """Applies the modifications described by a residue graft to a PDB residue belonging to the binder
    :param binder_residue: the affected PDB binder residue
    :param residue_graft: describes the modification to apply
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    """
    crep = residue_graft.pocket_cluster.get_most_representative_member(distance_factor, orient_factor, secor_factor)
    crep_icoor = get_icoor_from_aresidue(AbstractResidue(crep.binder_atoms))
    binres_icoor = get_icoor_from_aresidue(residue_graft.binder_res)
    binder_residue.resname = crep.binder_restype
    terminal = 'OXT' in binder_residue
    binder_residue.child_list = []
    binder_residue.child_dict = {}
    for atom_name, atom_position in crep.binder_atoms.items():
        if (crep.binder_restype == 'GLY' and atom_name == 'CB') or (not terminal and atom_name == 'OXT'):
            continue
        rel_position = crep_icoor.external_to_internal(atom_position, False)
        new_position = binres_icoor.internal_to_external(rel_position, False)
        atom = Atom(atom_name, new_position, 0.0, 1.0, ' ', atom_name, 0, atom_name[0])
        binder_residue.add(atom)


def apply_pocket_graft(ligand_residue: Residue, pocket_graft: PocketGraft,
                       distance_factor: float, orient_factor: float, secor_factor: float) -> None:
    """Updates the ligand coordinates based on the information obtained from a pocket graft.
    :param ligand_residue: the affected PDB ligand residue
    :param pocket_graft: this graft is made effective
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    """
    if len(pocket_graft.residue_grafts) == 0:
        return
    crep = pocket_graft.residue_grafts[0].pocket_cluster.get_most_representative_member(distance_factor, orient_factor,
                                                                                        secor_factor)
    icoor = get_icoor_from_aresidue(pocket_graft.ligand_res)
    ligand_residue.resname = crep.ligand_restype
    terminal = 'OXT' in ligand_residue
    ligand_residue.child_list = []
    ligand_residue.child_dict = {}
    for atom_name, atom_position in crep.ligand_atoms.items():
        if (crep.ligand_restype == 'GLY' and atom_name == 'CB') or (not terminal and atom_name == 'OXT'):
            continue
        if atom_name in ligand_residue:
            ligand_residue.detach_child(atom_name)
        new_position = icoor.internal_to_external(atom_position, False)
        atom = Atom(atom_name, new_position, 0.0, 1.0, ' ', atom_name, 0, atom_name[0])
        ligand_residue.add(atom)


def apply_binder_graft(scaffold: Structure, binder_graft: BinderGraft, distance_factor: float, orient_factor: float,
                       secor_factor: float, res_lock: List[Tuple[int, str]] = None) -> None:
    """Applies all pocket grafts and residue grafts contained in this binder graft to a given scaffold structure
    :param scaffold: the scaffold PDB structure that will be modified in-place.
    :param binder_graft: the binder graft to apply
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    :param res_lock: Locks binder residues and prevents them beeing mutated; List of Tuples containing residue id, type
    """

    def apply_residue_graft_in_presence_of_manual_mut():
        # If res_lock is defined, the scaffold defines some predefined mutations. Thus we have to check, if we would
        # destroy mutations with our residue graft. Plus mutations to same residue type is allowed, because predefined
        # mutations usually give no side chain information which we can add here.
        if (rg.binder_res.residue_id, newres) in res_lock or \
                rg.binder_res.residue_id not in (lock[0] for lock in res_lock):
            apply_residue_graft(res, rg, distance_factor, orient_factor, secor_factor)

    def apply_residue_graft_in_absence_of_manual_mut():
        # If res_lock is not defined, we did not define manual mutations to the scaffold. Thus, we don't have to check
        # the res_lock and we won't allow mutations to the original types.
        if newres != rg.binder_res.original_residue_type:
            print(f'{chain.get_id()}: {rg.binder_res.original_residue_type}{res.get_id()[1]}{newres}')
            apply_residue_graft(res, rg, distance_factor, orient_factor, secor_factor)

    if res_lock is None:
        apply_residue_graft_ = apply_residue_graft_in_absence_of_manual_mut
    else:
        apply_residue_graft_ = apply_residue_graft_in_presence_of_manual_mut

    for pg in binder_graft.pocket_grafts:
        for rg in pg.residue_grafts:
            for chain in scaffold.get_chains():
                if chain.get_id() == rg.binder_res.chain:
                    for res in chain.get_residues():
                        if res.get_id()[1] == rg.binder_res.residue_id:
                            newres = rg.pocket_cluster.get_most_representative_member(distance_factor, orient_factor,
                                                                                      secor_factor).binder_restype
                            apply_residue_graft_()
        for chain in scaffold.get_chains():
            if chain.get_id() == pg.ligand_res.chain:
                for res in chain.get_residues():
                    if res.get_id()[1] == pg.ligand_res.residue_id:
                        apply_pocket_graft(res, pg, distance_factor, orient_factor, secor_factor)


def simple_graft(scaffold: Structure, descriptor: ComplexDescriptor, pockets: Dict[str, List[Pocket]],
                 penalty_threshold: float, distance_factor: float, orient_factor: float, secor_factor: float,
                 select: Callable[[List[BinderGraft]], BinderGraft] = select_best_pocket_graft,
                 loc_observer: Callable[[int, int], None] = None,
                 glob_observer: Callable[[int, int], None] = None) -> Tuple[Union[BinderGraft, None], Structure]:
    """This function implements a simple algorithm for mutating the residues of a binder protein to match a ligand
    with given sequence. To find these mutations, information from the pocket miner is exploited.
    :param scaffold: the scaffold PDB structure that is modified in-place
    :param descriptor: describes the ligand residues as well as the wildcards of the binder chain
    :param pockets: result of the pocket miner where the algorithm extracts its information from
    :param penalty_threshold: see get_optimal_pocket_graft
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    :param select: the selection strategy to apply for choosing a pocket graft
    :param loc_observer: callback to react to the progress of the local optimization stage
    :param glob_observer: callback to react to the progress of the local optimization stage
    :return: the binder graft representing the best solution
    """
    design = copy.deepcopy(scaffold)
    acomplex = generate_complex(None, design, descriptor)

    opt_graft = get_optimal_graft_for_acomplex(acomplex=acomplex, pockets=pockets,
                                               penalty_threshold=penalty_threshold, distance_factor=distance_factor,
                                               orient_factor=orient_factor, secor_factor=secor_factor,
                                               select=select, loc_observer=loc_observer, glob_observer=glob_observer)
    if opt_graft is not None:
        apply_binder_graft(design, opt_graft, distance_factor, orient_factor, secor_factor)
    return opt_graft, design


def multi_graft(scaffold: Structure, descriptor: ComplexDescriptor, pockets: Dict[str, List[Pocket]],
                penalty_threshold: float, distance_factor: float, orient_factor: float, secor_factor: float,
                n_solutions: int = 20, observer: Callable[[int], None] = None) \
        -> List[Tuple[float, BinderGraft, Structure]]:
    """A variant of the simple graft algorithm; the top n pocket grafts are globally optimized and the different
    corresponding binder grafts are returned
    :param scaffold: the scaffold PDB structure that is modified in-place
    :param descriptor: describes the ligand residues as well as the wildcards of the binder chain
    :param pockets: result of the pocket miner where the algorithm extracts its information from
    :param penalty_threshold: see get_optimal_pocket_graft
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    :param n_solutions: maximum number of solutions to generate
    :param observer: callback to react to the index of the solution currently processed
    :return: a list of max size n_solutions, which contains the penalty, the binder graft, and the mutated structure
    belonging to every solution
    """

    def _select_ith_pocket_graft(candidates: List[BinderGraft]) -> BinderGraft:
        try:
            return candidates[i]
        except IndexError:
            pass

    result: List[Tuple[float, BinderGraft, Structure]] = []
    for i in range(n_solutions):
        binder_graft, structure = simple_graft(scaffold, descriptor, pockets, penalty_threshold, distance_factor,
                                               orient_factor, secor_factor, select=_select_ith_pocket_graft)
        if observer is not None:
            observer(i)
        if binder_graft is None:
            break
        result.append((binder_graft.penalty(), binder_graft, structure))
    return result


def get_mutations_single_p_graft(pocket: Pocket, scaffold: Structure, descriptor: ComplexDescriptor, res_id: int = 0,
                                 penalty_threshold: float = 32.0, distance_factor: float = 2.0,
                                 orient_factor: float = 1.0,
                                 secor_factor: float = 1.0, res_lock: List[Tuple[int, str]] = None) \
        -> Tuple[PocketGraft, Structure]:
    acomplex = generate_complex(None, scaffold, descriptor)
    pocket_grafts: List[PocketGraft] = []

    # Take ligand residue with correct residue id
    ligres = None
    for lres_i, lres in enumerate(acomplex.ligand.residues):
        if lres.residue_id == res_id:
            ligres = acomplex.ligand.residues[lres_i]
            break
    assert ligres is not None

    icoor = get_icoor_from_aresidue(ligres)
    wc_res = [br.transform_coor(icoor, True) for br in acomplex.binder.residues]
    pocket_graft = get_optimal_pocket_graft(ligres, wc_res, pocket, penalty_threshold, distance_factor,
                                            orient_factor, secor_factor)
    for i, rg in enumerate(pocket_graft.residue_grafts):
        pocket_graft.residue_grafts[i].binder_res = rg.binder_res.transform_coor(icoor, False)

    pocket_grafts.append(pocket_graft)
    binder_graft = BinderGraft(pocket_grafts)  # No need to find best graft, because we only have one
    apply_binder_graft(scaffold, binder_graft, distance_factor, orient_factor, secor_factor, res_lock=res_lock)
    return pocket_graft, scaffold


def prepare_scaffold_and_descriptor_for_grafting(scaffold_name: str,
                                                 ligand_residues: Union[Tuple[str, str],
                                                                        Dict[str, Dict[Union[str, int], str]]],
                                                 design_positions: Dict[str, List[Union[str, int]]],
                                                 pocket_restype: str = None):
    """
    Prepares the scaffold structure (Bio.PDB) from path string and the ComplexDescriptor from definitions
    of ligand residues and binder residues for subsequent grafting.

    :param scaffold_name: the path to a pdb structure including the chain ids and res ids given in ligand_residue
    and design_positions.
    :param ligand_residues:
    Either: A tuple of chain id, residue id (str or int) and desired residue type (needs a pocket_restype definition!)
    OR: A dict of chain id matched with a dict of residue id (str or int) and desired residue type.
    :param design_positions: A dictionary with the keys being the chain ids and the values being a list of all residue
    ids that should be allowed to be mutated while grafting.
    :param pocket_restype: The residue type the ligand should be mutated to. Only used if ligand_residues is a tuple.
    :return: scaffold as Bio.PDB structure and ComplexDescriptor
    """

    scaffold = PDBParser(QUIET=True).get_structure(scaffold_name,
                                                   scaffold_name)
    if isinstance(ligand_residues, dict):
        ligand_info = []
        for ligand_chain, ligand_chain_residues in ligand_residues.items():
            for ligand_res_id, ligand_restype in ligand_chain_residues.items():
                ligand_info.append((ligand_chain, str(ligand_res_id), ligand_restype))
    else:
        ligand_chain, ligand_res_id = ligand_residues
        ligand_info = [(ligand_chain, str(ligand_res_id), pocket_restype)]
        assert pocket_restype is not None
    binder_info = []
    for chain, positions in design_positions.items():
        binder_info += list((chain, str(pos)) for pos in positions)
    coupling_info = []
    descriptor = ComplexDescriptor(ligand_info, binder_info, coupling_info)

    return scaffold, descriptor


def graft_pocket_onto_scaffold(pocket: Pocket, scaffold_name: str, ligand_residue: Tuple[str, Union[str, int]],
                               design_positions: Dict[str, List[Union[str, int]]], output_name=None) -> Structure:
    """
    Will graft the residues of pocket onto the pdb structure found at scaffold_name.
    The grafting is applied based on the ligand residue position and the design_positions (where mutations are
    permitted).

    :param pocket: a Pocket object containing at least one member per cluster
    :param scaffold_name: the path to a pdb structure including the chain ids and res ids given in ligand_residue
    and design_positions.
    :param ligand_residue: A tuple of chain id, residue id (str or int) and desired residue type.
    :param design_positions: A dictionary with the keys being the chain ids and the values being a list of all residue
    ids that should be allowed to be mutated while grafting.
    :param output_name: If given this will be the path of the newly designed pdb structure
    :return: The newly designed Bio.PDB structure
    """
    scaffold, descriptor = prepare_scaffold_and_descriptor_for_grafting(scaffold_name=scaffold_name,
                                                                        ligand_residues=ligand_residue,
                                                                        design_positions=design_positions,
                                                                        pocket_restype=pocket.restype)

    graft, design = get_mutations_single_p_graft(pocket, scaffold, descriptor, res_id=int(ligand_residue[1]))
    if output_name is not None:
        io = PDBIO()
        io.set_structure(design)
        io.save(output_name)

    return design


def graft_best_pocket_onto_scaffold(pockets: Dict[str, List[Pocket]], scaffold_name: str,
                                    ligand_residues: Dict[str, Dict[Union[str, int], str]],
                                    design_positions: Dict[str, List[Union[str, int]]], output_name=None,
                                    penalty_threshold: float = 32.0, distance_factor: float = 2.0,
                                    orient_factor: float = 1.0, secor_factor: float = 1.0, ):
    """
    Will graft the residues of the best matching pocket in pockets onto the pdb structure found at scaffold_name.
    The grafting is applied based on the ligand residue positions and the design_positions (where mutations are
    permitted) and based on the defined ligand residue types.

    :param pockets: a collection of Pockets (result of the pocket miner) containing at least one member per cluster
    and at least one pocket where the restype is the ligand_restype.
    :param scaffold_name: the path to a pdb structure including the the chain ids and res ids given in ligand_residue
    and design_positions.
    :param ligand_residues: A dict of chain id matched with a dict of residue id (str or int) and desired residue type.
    :param design_positions: A dictionary with the keys being the chain ids and the values being a list of all residue
    ids that should be allowed to be mutated while grafting.
    :param output_name: If given this will be the path of the newly designed pdb structure
    :param penalty_threshold: see get_optimal_pocket_graft
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    :return: A tuple containing the newly designed Bio.PDB structure and the BinderGraft applied
    """
    scaffold, descriptor = prepare_scaffold_and_descriptor_for_grafting(scaffold_name=scaffold_name,
                                                                        ligand_residues=ligand_residues,
                                                                        design_positions=design_positions)

    graft, design = simple_graft(scaffold=scaffold, descriptor=descriptor, pockets=pockets,
                                 penalty_threshold=penalty_threshold, distance_factor=distance_factor,
                                 orient_factor=orient_factor, secor_factor=secor_factor)
    if output_name is not None:
        io = PDBIO()
        io.set_structure(design)
        io.save(output_name)

    return design, graft


def multi_graft_best_pocket_onto_scaffold(pockets: Dict[str, List[Pocket]], scaffold_name: str,
                                          ligand_residues: Dict[str, Dict[Union[str, int], str]],
                                          design_positions: Dict[str, List[Union[str, int]]], output_name=None,
                                          penalty_threshold: float = 32.0, distance_factor: float = 2.0,
                                          orient_factor: float = 1.0, secor_factor: float = 1.0,
                                          n_solutions: int = 10):
    """
    Does the same thing as graft_best_pocket_onto_scaffold, but returns the best n_solutions solutions of grafting:

    Will graft the residues of the best matching pocket in pockets onto the pdb structure found at scaffold_name.
    The grafting is applied based on the ligand residue positions and the design_positions (where mutations are
    permitted) and based on the defined ligand residue types.

    :param pockets: a collection of Pockets (result of the pocket miner) containing at least one member per cluster
    and at least one pocket where the restype is the ligand_restype.
    :param scaffold_name: the path to a pdb structure including the the chain ids and res ids given in ligand_residue
    and design_positions.
    :param ligand_residues: A dict of chain id matched with a dict of residue id (str or int) and desired residue type.
    :param design_positions: A dictionary with the keys being the chain ids and the values being a list of all residue
    ids that should be allowed to be mutated while grafting.
    :param output_name: If given this will be the path of the newly designed pdb structures with added numbers for
    the n_th solution.
    :param penalty_threshold: see get_optimal_pocket_graft
    :param distance_factor: see calculate_penalty
    :param orient_factor: see calculate_penalty
    :param secor_factor: see calculate_penalty
    :param n_solutions: maximum number of solutions to generate
    :return: A list containing: penalty, BinderGraft and designed Bio.PDB structures (sorted best to worst).
    """
    scaffold, descriptor = prepare_scaffold_and_descriptor_for_grafting(scaffold_name=scaffold_name,
                                                                        ligand_residues=ligand_residues,
                                                                        design_positions=design_positions)

    result = multi_graft(scaffold=scaffold, descriptor=descriptor, pockets=pockets,
                         penalty_threshold=penalty_threshold, distance_factor=distance_factor,
                         orient_factor=orient_factor, secor_factor=secor_factor, n_solutions=n_solutions)

    if output_name is not None:
        if not output_name.endswith(".pdb"):
            output_name = output_name + ".pdb"
        for number, result_n in enumerate(result):
            design = result_n[2]
            io = PDBIO()
            io.set_structure(design)
            io.save(output_name.replace(".pdb", str(number + 1) + ".pdb"))

    return result
