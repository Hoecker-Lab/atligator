"""This module provides data strucures and algoriths for mining frequent host pockets binding to specific ligand
residue types based on atlas data. It is based on association rule mining and frequent itemset mining algorithms.

:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-04-20
"""
import json
import os
from concurrent.futures import ThreadPoolExecutor
from copy import copy
from math import sqrt
from random import sample, choice
from typing import List, Dict, Tuple, Set, Callable, TextIO, Union, Any

from Bio.PDB import Vector
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.SeqUtils import seq1
from scipy.optimize import linear_sum_assignment

from atligator.atlas import Atlas, AtlasDatapoint, AtlasEncoder
from atligator.pdb_util import canonical_amino_acids, find_scop_i_from_pdbs


class PocketItemset:
    """A pocket itemset is a collection of binder residues found to be interacting with a specific ligand residue."""
    items: Dict[str, int]

    def __init__(self, items: Dict[str, int] or str):
        """Creates a new instance.
        :param items: relates the types of items (i.e. their binder residue type) with the number of occurrences in
        a pocket. 'XXX' is a shorthand for the dict {'XXX': 1}
        """
        if type(items) is str:
            self.items = {items: 1}
        else:
            self.items = items

    def pocket_id(self) -> str:
        """
        :return: a sequence of one-letter codes of items part of the atlas. Multiple containment is expressed by
        duplication of letters. Example: 'KRRQNSS'
        """
        result = ''
        for k, v in self.items.items():
            for _ in range(v):
                result += seq1(k)
        return result

    def __contains__(self, other: 'PocketItemset') -> bool:
        """
        :param other:
        :return: whether this pocket itemset contains all items (and as least as many of duplicates) as other.
        """
        for itype, icount in other.items.items():
            if itype not in self.items or self.items[itype] < icount:
                return False
        return True

    def __add__(self, other: 'PocketItemset') -> 'PocketItemset':
        """
        :param other:
        :return: pocket itemset that represents the union of this and other.
        """
        u_items: Dict[str, int] = {}
        for itype in list(self.items.keys()) + list(other.items.keys()):
            icount_self = 0 if itype not in self.items else self.items[itype]
            icount_other = 0 if itype not in other.items else other.items[itype]
            u_items[itype] = icount_self + icount_other
        return PocketItemset(u_items)

    def __eq__(self, other) -> bool:
        return self.items == other.items

    def __hash__(self) -> int:
        return hash(frozenset(self.items.items()))

    def __len__(self) -> int:
        return sum(v for v in self.items.values())

    def __str__(self) -> str:
        return str(self.items)

    def __repr__(self) -> str:
        return str(self)


class PocketItemDataset:
    """A pocket item dataset organizes a collection of pocket itemsets and allows to make statistical statements
    (e.g., support, confidence). Furthermore, it enables the identification of frequent itemsets."""

    def __init__(self, itemsets: List[PocketItemset]):
        """Creates a new instance of the pocket dataset.
        :param itemsets: the initial collection of pocket itemsets
        """
        self.itemsets = itemsets
        self._support_cache: Dict[PocketItemset, float] = {}

    def implied_itemsets(self, itemset: PocketItemset) -> List[PocketItemset]:
        """Retrieves all itemsets managed by this dataset which contain the items described by the parameter itemset.
        :param itemset:
        """
        return list([i for i in self.itemsets if itemset in i])

    def support(self, itemset: PocketItemset) -> float:
        """The support of an itemset is defined as the ratio between the number of those itemsets that imply it and
        the number of itemsets managed in total.
        For performance reason, returned values are cached in the _support_cache, which is transparent to the caller.
        :param itemset:
        :return: the support of itemset. Ranges between 0 and 1.
        """
        if itemset in self._support_cache:
            return self._support_cache[itemset]
        else:
            lis = len(self.itemsets)
            support = 0.0 if lis == 0 else float(len(self.implied_itemsets(itemset))) / float(lis)
            self._support_cache[itemset] = support
            return support

    def confidence(self, premise: PocketItemset, conclusion: PocketItemset):
        """The confidence of an association rule X->Y between two itemsets is defined as the support of the union of X
        and Y, divided by the support of X.
        :param premise: left hand side of the association rule, i.e., X
        :param conclusion: right hand side, Y
        :return: the confidence of X->Y. Ranges between 0 and 1.
        """
        supp_pre = self.support(premise)
        return 0.0 if supp_pre == 0.0 else self.support(premise + conclusion) / supp_pre

    def find_associations(self, itemset: PocketItemset, confidence_threshold: float) -> Dict[str, float]:
        """Extends the given itemset by every possible item, taken from the list of 20 canonical amino acids. Relates
        the added item to the confidence related to the addition of the new item
        :param itemset: the base pocket itemset to find new association rules for
        :param confidence_threshold: associations showing a confidence below this value will not be returned
        :return: dict that maps potential residue types to the confidence predicted
        """
        associations: Dict[str, float] = {}
        # TODO alternative aas -> atm: associations related to binder residues, thus, not needed
        for restype in canonical_amino_acids:
            confidence = self.confidence(itemset, PocketItemset(restype))
            if confidence > confidence_threshold:
                associations[restype] = confidence
        return associations

    def find_frequent_itemsets(self, confidence_threshold: float, support_threshold: float,
                               min_itemset_cardinality: int,
                               cardinality_base: float = 2.0) -> List[Tuple[float, PocketItemset]]:
        """Uses a breadth-first search to internally create a frequent itemset lattice based on repeated applications
        of the find_associations function. Returns frequent itemsets associated with their support.
        :param confidence_threshold: association rules that fall below this threshold are not applied
        :param support_threshold: itemsets below this value are not considered as frequent
        :param min_itemset_cardinality: itemsets whose cardinality is equal or greater than this value are included in
        the result
        :param cardinality_base: used for fitness calculation of each itemset in 'support * sqrt(base)^|itemset|'.
        Higher values result in pockets with more items.
        :return: a list of tuples relating support values to corresponding itemsets mined from this dataset
        """
        result: List[Tuple[float, PocketItemset]] = []
        tabu: Set[PocketItemset] = set()
        previous: List[PocketItemset] = [PocketItemset({})]
        while len(previous) > 0:
            current = []
            for iset in previous:
                for new_item in self.find_associations(iset, confidence_threshold).keys():
                    new_itemset = iset + PocketItemset(new_item)
                    if new_itemset not in tabu:
                        tabu.add(new_itemset)
                        new_supp = self.support(new_itemset)
                        if new_supp > support_threshold:
                            current.append(new_itemset)
                            if len(new_itemset) >= min_itemset_cardinality:
                                result.append((new_supp, new_itemset))
            previous = current
        return sorted(result, key=lambda e: e[0] * float(sqrt(cardinality_base) ** float(len(e[1]))), reverse=True)


def extract_pocket_item_datasets(atlas: Atlas) -> Dict[str, PocketItemDataset]:
    """Evaluates an atlas and returns up to 20 pocket item datasets, one for each residue type occurring in the atlas.
    :param atlas: the input for the analysis
    :return: dict from ligand residue type to their individual item dataset
    """
    item_datasets: Dict[str, PocketItemDataset] = {}
    for lig_origin, datapoint_group in atlas.group_by_ligand_origin().items():
        itemset = PocketItemset({})
        for datapoint in datapoint_group:
            itemset += PocketItemset(datapoint.binder_restype)
        lig_restype = datapoint_group[0].ligand_restype
        if lig_restype not in item_datasets:
            item_datasets[lig_restype] = PocketItemDataset([])
        item_datasets[lig_restype].itemsets.append(itemset)
    return item_datasets


def filter_atlas_by_itemset(atlas: Atlas, ligand_restype: str, itemset: PocketItemset,
                            include_unrelated_residues: bool = False) -> Atlas:
    """Restricts a given atlas to those ligand residues whose binding pocket is represented by a given itemset. In
    this way, an atlas for a specific binding pocket can be obtained automatically.
    :param atlas: the atlas to filter
    :param ligand_restype: ligand residue type to assume for filtering. None if this should not be used as a criterion
    (in this case, the resulting atlas would still represent several ligand residue types)
    :param itemset: datapoint groups that match this itemset pass the filter
    :param include_unrelated_residues: if true, the resulting atlas also includes residue types that are not part of
    the actual pocket. This introduces more noise, but can also indicate interactions overlooked by itemset mining
    :return: the filtered atlas
    """
    filtered_dps: List[AtlasDatapoint] = []
    for _, datapoints in atlas.group_by_ligand_origin().items():
        g_lig_restype = datapoints[0].ligand_restype
        if ligand_restype is not None and g_lig_restype != ligand_restype:
            continue
        g_itemset = PocketItemset({})
        for datapoint in datapoints:
            g_itemset += PocketItemset(datapoint.binder_restype)
        if itemset in g_itemset:
            for datapoint in datapoints:
                if include_unrelated_residues or datapoint.binder_restype in itemset.items.keys():
                    filtered_dps.append(datapoint)
    return Atlas(filtered_dps)


def distance(res1: AtlasDatapoint, res2: AtlasDatapoint,
             distance_factor: float, orient_factor: float, secor_factor: float) -> float:
    """Distance function used for clustering of residues. Weighted sum of the distances between C_alpha vectors, and
     primary and secondary orientation angles.
    :param res1: first operand
    :param res2: second operand
    :param distance_factor: weight for the distance between c_alpha positions
    :param orient_factor: weight for the difference between primary orientation
    :param secor_factor: weight for the difference between secondary orientation
    :return: the distance calculated (a nonnegative scalar)
    """
    calpha_distance: float = abs((res1.binder_calpha() - res2.binder_calpha()).norm())
    orient_angle: float = abs(res1.binder_orientation().angle(res2.binder_orientation()))
    sec_orient_angle: float = abs(res1.binder_secondary_orientation().angle(res2.binder_secondary_orientation()))
    return distance_factor * calpha_distance + orient_factor * orient_angle + secor_factor * sec_orient_angle


class PocketCluster:
    """A cluster of atlas residues that represents a part of a pocket."""

    @classmethod
    def from_json(cls, cluster: Dict[str, Any]):
        """
        Return a PocketCluster from a json representation of one.
        :param cluster: json representation of a PocketCluster or open file handle containing it
        :return: PocketCluster based on json
        """
        centroid = AtlasDatapoint.from_json(cluster['centroid'])
        new_cluster = PocketCluster(ligand_restype=cluster['ligand_restype'],
                                    binder_restype=cluster['binder_restype'],
                                    initial_centroid=centroid)
        new_cluster.variance = cluster['variance']
        new_cluster.has_changed = cluster['has_changed']
        new_cluster.members = [AtlasDatapoint.from_json(mem) for mem in cluster['members']]
        new_cluster.future_members = [AtlasDatapoint.from_json(fmem) for fmem in cluster['future_members']]
        return new_cluster

    def __init__(self, ligand_restype: str, binder_restype: str, initial_centroid: AtlasDatapoint):
        """ Creates a new cluster.
        :param ligand_restype: same as for all clustered atlas residues
        :param binder_restype: same as for all clustered atlas residues
        :param initial_centroid: the initial centroid of the cluster. Not counted as a member.
        """
        self.ligand_restype = ligand_restype
        self.binder_restype = binder_restype
        self.members: List[AtlasDatapoint] = []
        self.centroid = copy(initial_centroid)
        self.variance = 0.0
        self.future_members: List[AtlasDatapoint] = []
        self.has_changed = True

    def __len__(self):
        return len(self.members)

    def get_centroid(self) -> AtlasDatapoint:
        """Calculates a new centroid based on the members of this clusters. It is obtained by calculating the
        arithmetic means of the components of all c_alpha, c_beta, and c_o vectors.
        :return: an artificial atlas residue that represents the centroid of this cluster
        """
        cent = AtlasDatapoint(self.ligand_restype, self.binder_restype, None, None, {}, {})
        cl = 1.0 / float(len(self.future_members))
        for ares in self.future_members:
            for lig_atom, lig_pos in ares.ligand_atoms.items():
                if lig_atom not in cent.ligand_atoms:
                    cent.ligand_atoms[lig_atom] = Vector(0.0, 0.0, 0.0)
                cent.ligand_atoms[lig_atom] += ares.ligand_atoms[lig_atom] ** cl
            for bin_atom, bin_pos in ares.binder_atoms.items():
                if bin_atom not in cent.binder_atoms:
                    cent.binder_atoms[bin_atom] = Vector(0.0, 0.0, 0.0)
                cent.binder_atoms[bin_atom] += ares.binder_atoms[bin_atom] ** cl
        return cent

    def update(self, distance_factor: float, orient_factor: float, secor_factor: float, keep_members: bool = False) \
            -> None:
        """To be called after every clustering iteration. Updates the centroid, the cached variance, saves the list
        of members in order to detect whether the contents have changed, and clears the current member list.
        :param distance_factor: weight for the distance between c_alpha positions in distance function
        :param orient_factor: weight for the difference between primary orientation in distance function
        :param secor_factor: weight for the difference between secondary orientation in distance function
        """
        if len(self.future_members) > 0:
            self.centroid = self.get_centroid()
        self.variance = self.get_variance(None, distance_factor, orient_factor, secor_factor)
        self.has_changed = frozenset(self.members) != frozenset(self.future_members)
        if keep_members:
            self.members = self.members + self.future_members[:]
        else:
            self.members = self.future_members[:]
        self.future_members.clear()

    def get_variance(self, additional_member: AtlasDatapoint or None,
                     distance_factor: float, orient_factor: float, secor_factor: float) -> float:
        """Obtains the variance of this cluster, which is obtained by the sum of squares of distances between all
        members and the centroid, divided by the number of members.
        :param additional_member: optionally, a hypothetical member can be provided here. It will be treated as if it
        were actually part of the cluster. Otherwise, fill this parameter with None.
        :param distance_factor: weight for the distance between c_alpha positions in distance function
        :param orient_factor: weight for the difference between primary orientation in distance function
        :param secor_factor: weight for the difference between secondary orientation in distance function
        :return: the variance of this cluster, a nonnegative scalar. Square of the standard deviation
        """
        if len(self.future_members) == 0 and additional_member is None:
            return 0.0
        dist = 0.0
        for ares in self.future_members + ([] if additional_member is None else [additional_member]):
            dist += (distance(self.centroid, ares, distance_factor, orient_factor, secor_factor) ** 2.0)
        return dist / float(len(self.future_members) + (0 if additional_member is None else 1))

    def get_variance_increase(self, candidate: AtlasDatapoint,
                              distance_factor: float, orient_factor: float, secor_factor: float) -> float:
        """Calculates the increase in variance that would be obtained by adding a given candidate to the list of members
        of this cluster.
        :param candidate: virtual additional member of the cluster (will not be added)
        :param distance_factor: weight for the distance between c_alpha positions in distance function
        :param orient_factor: weight for the difference between primary orientation in distance function
        :param secor_factor: weight for the difference between secondary orientation in distance function
        :return: variance including the candidate minus variance without the candidate. May be negative
        """
        return self.get_variance(candidate, distance_factor, orient_factor, secor_factor) - self.get_variance(
            None, distance_factor, orient_factor, secor_factor)

    def get_least_representative_member(self, distance_factor: float, orient_factor: float, secor_factor: float) \
            -> AtlasDatapoint:
        """
        :return: the member of the cluster that has the largest distance to the centroid.
        :param distance_factor: weight for the distance between c_alpha positions in distance function
        :param orient_factor: weight for the difference between primary orientation in distance function
        :param secor_factor: weight for the difference between secondary orientation in distance function
        """
        max_dist = 0.0
        lrm = None
        for m in self.members:
            dist = distance(self.centroid, m, distance_factor, orient_factor, secor_factor)
            if dist > max_dist:
                max_dist = dist
                lrm = m
        return lrm

    def get_most_representative_member(self, distance_factor: float, orient_factor: float, secor_factor: float) \
            -> AtlasDatapoint:
        """
        :return: the member of the cluster that has the shortest distance to the centroid.
        :param distance_factor: weight for the distance between c_alpha positions in distance function
        :param orient_factor: weight for the difference between primary orientation in distance function
        :param secor_factor: weight for the difference between secondary orientation in distance function
        """
        min_dist = float('inf')
        mrm = None
        for m in self.members:
            dist = distance(self.centroid, m, distance_factor, orient_factor, secor_factor)
            if dist < min_dist:
                min_dist = dist
                mrm = m
        return mrm

    def get_dist_to_centroid(self, other: AtlasDatapoint, distance_factor: float,
                             orient_factor: float, secor_factor: float) -> float:
        return distance(self.centroid, other, distance_factor, orient_factor, secor_factor)


class Pocket:
    """A pocket is a collection of pocket clusters that have emerged from a given itemset."""

    @classmethod
    def from_json(cls, pocket):
        """
        Return a Pocket from a json representation of one.
        :param pocket: json representation of a Pocket
        :return: Pocket based on json
        """
        pocket_itemset = PocketItemset(items=pocket['itemset'])
        pocket_clusters = []
        for cluster in pocket['clusters']:
            pocket_clusters.append(PocketCluster.from_json(cluster))
        return Pocket(restype=pocket['restype'],
                      support=pocket['support'],
                      itemset=pocket_itemset,
                      clusters=pocket_clusters)

    def to_json(self):
        """
        Return a json representation of this Pocket.
        :return: json string based on Pocket
        """
        return json.dumps(self, cls=PocketsEncoder)

    def __init__(self, restype: str, itemset: PocketItemset, support: float, clusters: List[PocketCluster]):
        """
        :param restype: the ligand residue type represented by this pocket
        :param itemset: the itemset containing the binder residue types part of the binding pocket
        :param support: absolute support of the pocket's itemset in the dataset
        :param clusters: contains the positions of clustered atlas residues for every element of the itemset at
        corresponding positions.
        """
        self.restype = restype
        self.itemset = itemset
        self.support = support
        self.clusters = clusters

    def __repr__(self):
        return f"<Pocket {self.restype}: {self.itemset}>"

    def __len__(self) -> int:
        return self.n_ligand_origins()

    def pocket_id(self) -> str:
        """
        :return: an identifier for this pocket, based on its underlying dataset
        """
        return self.itemset.pocket_id()

    def n_ligand_origins(self) -> int:
        """
        :return: the number of distinct ligand origins the members of this cluster have emerged from
        """
        lig_ors: Set[str] = set()
        for c in self.clusters:
            for m in c.members:
                lig_ors.add(m.ligand_origin)
        return len(lig_ors)

    def origins(self, scop_db: str = None,
                stringify: bool = True, scop_only: bool = True) -> Dict or str or None:
        """
        Specfies origins of data points in this pocket (in scope classifiers)
        :param scop_db: Scope database reference
        :param stringify: If true the output will be a string
        :param scop_only: If true the output will be a dictionary with number scope occurances (or its string)
        :return: Either a dict of scope identifiers with number of occurances or
        a dict with pdb names and their correlated scope identifiers
        """

        def get_scop_i(pdb_scop_dict, sort_d: bool = True) -> Dict or List:
            """
            Returns number of occurances of this pocket in a specific SCOPe family
            :param pdb_scop_dict: Dictionary to find pdb - scop connections from
            :param sort_d: If true the output will be a list of keys and values of the output dictionary
            sorted by number of occurances
            :return: Dictionary with scope identifiers and the number of occurances of this pocket
            """
            scop_i_n: Dict[str, int] = {"None": 0}
            for l in pdb_scop_dict.values():
                if l is None:
                    scop_i_n["None"] += 1
                    continue
                for i in l:
                    if i not in scop_i_n.keys():
                        scop_i_n[i] = 1
                    else:
                        scop_i_n[i] += 1

            if not sort_d:
                return scop_i_n
            return [(k, scop_i_n[k]) for k in sorted(scop_i_n, key=scop_i_n.__getitem__, reverse=True)]

        if os.path.isfile(scop_db):
            pdb_set: Set[str] = set()
            for c in self.clusters:
                for m in c.members:
                    pdb_set.add(m.ligand_origin)
            pdbs: List[str] = []
            for pdb in pdb_set:
                pdbs.append(pdb[:4])
            d = find_scop_i_from_pdbs(pdbs, scop_db)
            if scop_only:
                return get_scop_i(d)
            return d
        return None

    def to_atlas(self) -> Atlas:
        """Converts a pocket into an atlas based on the members of the clusters.
        :return: an atlas that represents the pocket
        """
        a_residues: List[AtlasDatapoint] = []
        for c in self.clusters:
            a_residues.extend(c.members)
        return Atlas(a_residues)

    def get_most_representative_group(self, distance_factor: float, orient_factor: float, secor_factor: float):
        """ Selects the most representative group of cluster members. A group consists of a ligand residue and
        binder residues. The number of binder residues corresponds to the number of clusters and thus to the cardinality
        of this pocket.
        :param distance_factor: weight for the distance between c_alpha positions in distance function
        :param orient_factor: weight for the difference between primary orientation in distance function
        :param secor_factor: weight for the difference between secondary orientation in distance function
        :return: an atlas that represents the clustered pocket
        """
        min_len = 2 ** 10000
        min_c = None
        # TODO only applicable if all clusters must include at least one..
        for c_i, c in enumerate(self.clusters):  # Find smallest cluster first
            if len(c.members) < min_len:
                min_c = c_i
                min_len = len(c.members)
        suited_origins: Dict[str, List] = {}
        for member in self.clusters[min_c].members:  # Find the ligand origins that can be found in all clusters
            if all(member.ligand_origin in (m.ligand_origin for m in c.members) for c in self.clusters):
                suited_origins[member.ligand_origin] = [None for _ in range(len(self.clusters))]
        if not len(suited_origins):  # If no suited origin is found, no group can be assigned
            print(self.pocket_id())
            return None, None

        origin_sums = {}
        for so in suited_origins.keys():  # For each suited origins find the best representatives in each cluster
            for c_i, c in enumerate(self.clusters):
                minimum = float('inf')
                for m in c.members:
                    if m.ligand_origin == so:
                        dist = c.get_dist_to_centroid(m, distance_factor, orient_factor, secor_factor)
                        if dist < minimum:
                            suited_origins[so][c_i] = [m.binder_origin, dist]  # Save the binder origin
            origin_sums[so] = sum(d[1] for d in suited_origins[so])
        best_origin = min(origin_sums, key=origin_sums.get)
        return suited_origins[best_origin], best_origin

    def to_clustered_atlas(self, distance_factor: float, orient_factor: float, secor_factor: float,
                           most_represent_group: bool = False) -> Atlas:
        """Converts a pocket into an atlas based on the centroids of the clusters.
        :param distance_factor: weight for the distance between c_alpha positions in distance function
        :param orient_factor: weight for the difference between primary orientation in distance function
        :param secor_factor: weight for the difference between secondary orientation in distance function
        :return: an atlas that represents the clustered pocket
        """
        ca_residues: List[AtlasDatapoint] = []
        if most_represent_group:  # Select for best grouped binder residues within this itemset
            members_of_best_lig_origin, best_ligand_origin = \
                self.get_most_representative_group(distance_factor, orient_factor, secor_factor)
            if best_ligand_origin is None:
                return Atlas([])
            for c_i, c in enumerate(self.clusters):
                for member in c.members:
                    if member.binder_origin == members_of_best_lig_origin[c_i][0] \
                            and member.ligand_origin == best_ligand_origin:
                        ca_residues.append(member)
                        break

        else:
            for bin_restype, count in self.itemset.items.items():
                r_clusters = list([c for c in self.clusters if c.binder_restype == bin_restype])
                for c in r_clusters:
                    ca_residues.append(c.get_most_representative_member(distance_factor, orient_factor, secor_factor))
        return Atlas(ca_residues)

    def to_pdb_structure(self) -> Structure:
        """
        :return: A PDB compliant structure that contains the datapoints of this pocket. By convention, chain A contains
        different conformations of the ligand residues of all datapoints, chains B and further contain residue
        coordinates for the conformations represented by individual members of the clusters.
        """
        model = Model(0)
        resid = 1
        ligand_chain = Chain('A')
        ligand_origins: Set[str] = set()
        for cluster in self.clusters:
            for member in cluster.members:
                if member.ligand_origin in ligand_origins:
                    continue
                ligres = member.ligand_to_pdb_residue(resid)
                ligand_chain.add(ligres)
                ligand_origins.add(member.ligand_origin)
                resid += 1
        model.add(ligand_chain)
        for i, cluster in enumerate(self.clusters):
            binder_chain_id = chr(ord('B') + i)
            binder_chain = Chain(binder_chain_id)
            for member in cluster.members:
                binres = member.binder_to_pdb_residue(resid)
                binder_chain.add(binres)
                resid += 1
            model.add(binder_chain)
        struc = Structure(self.pocket_id())
        struc.add(model)
        return struc

    def to_clustered_pdb_structure(self, distance_factor: float, orient_factor: float, secor_factor: float,
                                   most_represent_group: bool = True) -> Structure:
        """
        :param distance_factor: weight for the distance between c_alpha positions in distance function
        :param orient_factor: weight for the difference between primary orientation in distance function
        :param secor_factor: weight for the difference between secondary orientation in distance function
        :param most_represent_group: If True, the best group of residues (originating from one pdb structure!) are
        selected. Otherwise, the best individual residues are selected even if the originate from different structures
        and might not be senseful in combination (sterical hindering etc.).

        :return: A PDB compliant structure that contains the data of this pocket in a clustered way. By convention,
        chain A contains the best-matching (i.e. closest to centroid) conformation of the ligand residues of all
        datapoints, chains B and further contain best-matching residue coordinates for the conformations represented
        by individual members of the clusters. Thus, every chain posesses one residue.
        """
        model = Model(0)
        resid = 1
        ligand_chain = Chain('A')
        repr_ligand_origins: Set[str] = set()
        if most_represent_group:  # Select for best grouped binder residues within this itemset
            members_of_best_lig_origin, best_ligand_origin = \
                self.get_most_representative_group(distance_factor, orient_factor, secor_factor)
            if best_ligand_origin is None:
                return Structure(self.pocket_id())
            for c in self.clusters:
                for member in c.members:
                    if member.ligand_origin == best_ligand_origin:
                        ligres = member.ligand_to_pdb_residue(resid)
                        ligand_chain.add(ligres)
                        resid += 1
                        break
                break
            chain_counter = 0
            for c_i, c in enumerate(self.clusters):
                for member in c.members:
                    if member.binder_origin == members_of_best_lig_origin[c_i][0] \
                            and member.ligand_origin == best_ligand_origin:
                        binres = member.binder_to_pdb_residue(resid)
                        binder_chain_id = chr(ord('B') + chain_counter)
                        binder_chain = Chain(binder_chain_id)
                        binder_chain.add(binres)
                        resid += 1
                        chain_counter += 1
                        model.add(binder_chain)
                        break
            model.add(ligand_chain)
            struc = Structure(self.pocket_id())
            struc.add(model)
        else:  # find the best fitting individual residues within this itemset
            for cluster in self.clusters:
                cluster_repr = cluster.get_most_representative_member(distance_factor, orient_factor, secor_factor)
                if cluster_repr is None or cluster_repr.ligand_origin in repr_ligand_origins:
                    continue
                ligres = cluster_repr.ligand_to_pdb_residue(resid)
                ligand_chain.add(ligres)
                repr_ligand_origins.add(cluster_repr.ligand_origin)
                resid += 1
                break
            for i, cluster in enumerate(self.clusters):
                cluster_repr = cluster.get_most_representative_member(distance_factor, orient_factor, secor_factor)
                if cluster_repr is None:
                    continue
                binder_chain_id = chr(ord('B') + i)
                binder_chain = Chain(binder_chain_id)
                binres = cluster_repr.binder_to_pdb_residue(resid)
                binder_chain.add(binres)
                resid += 1
                model.add(binder_chain)
            model.add(ligand_chain)
            struc = Structure(self.pocket_id())
            struc.add(model)
        return struc


def pockets_to_json_file(pockets: Dict[str, List[Pocket]], fp):
    """
    Dumps the json represenation of Pockets to an open file handle
    :param pockets: Dict with Pockets to convert to json
    :param fp: open file handle to write json into
    :return: None
    """
    return json.dump(pockets, fp=fp, cls=PocketsEncoder)


def pockets_to_json(pockets: Dict[str, List[Pocket]]):
    """
    Returns the json represenation of Pockets
    :param pockets: Dict with Pockets to convert to json
    :return: json represenation of Pockets
    """
    return json.dumps(pockets, cls=PocketsEncoder)


def json_to_pockets(json_str_or_file: Union[str, TextIO]):
    """
    Return a Pockets dict from a json representation of one.
    :param json_str_or_file: json representation of Pockets or open file handle containing it
    :return: dictionary with Pockets
    """

    def convert(pockets: Dict[str, List[Dict]]):
        return_pockets = {}
        for ligand_aa, ligand_pockets in pockets.items():
            return_pockets[ligand_aa] = []
            for pocket in ligand_pockets:
                return_pockets[ligand_aa].append(Pocket.from_json(pocket))
        return return_pockets

    if isinstance(json_str_or_file, str):
        return convert(json.loads(json_str_or_file))
    else:
        return convert(json.load(json_str_or_file))


def clusterize_pocket_atlas(atlas: Atlas, itemset: PocketItemset, support: float, variance_threshold: float,
                            distance_factor: float, orient_factor: float, secor_factor: float,
                            one_per_cluster: bool = True, more_than_one_ok: bool = True) -> Pocket:
    """Takes an atlas, which already represents an individual binding pocket, and an itemset this atlas must be based
    on. Calculates a list of pocket clusters, which non-disjointly divide up the residues of the input atlas. For every
    residue type, at least as many pocket clusters as reflected by the itemset are obtained.
    The implementation is based on Lloyd's k-means algorithm. In addition, a post-processing step ensures that the
    obtained clusters' variance is below a given variance threshold by removing outliers iteratively.
    :param atlas: the pocket atlas to clusterize
    :param itemset: the itemset the pocket atlas is based on
    :param support: support of the itemset
    :param variance_threshold: in case the calculated clusters display a variance above this threshold, they will be
    divided up hierarchically in an iterative post-processing step.
    :param distance_factor: weight for the distance between c_alpha positions in distance function
    :param orient_factor: weight for the difference between primary orientation in distance function
    :param secor_factor: weight for the difference between secondary orientation in distance function
    :return: a list of pocket clusters representing a non-disjoint decomposition of the pocket atlas
    """

    def __assign_to_best_cluster(candidate_clusters: List[PocketCluster], datapoint: AtlasDatapoint) -> None:
        """Utility function that finds the most suitable among a set of candidate clusters, such that the variance of
        the cluster is increased to the lowest possible extent. The given datapoint is then appended to this cluster.
        :param candidate_clusters: the set of clusters to find suitable candidates within
        :param datapoint: the datapoint to decide cluster membership for
        """
        min_vi = float('inf')
        min_cluster: PocketCluster = None
        for c in candidate_clusters:
            if c.binder_restype == datapoint.binder_restype:
                vi = c.get_dist_to_centroid(datapoint, distance_factor, orient_factor, secor_factor)
                if vi < min_vi:
                    min_vi = vi
                    min_cluster = c
        min_cluster.future_members.append(datapoint)

    def __get_restypes(datapoints: List[AtlasDatapoint]):
        restypes = {}
        for dp_ in datapoints:
            try:
                restypes[dp_.binder_restype] += 1
            except KeyError:
                restypes[dp_.binder_restype] = 1
        return restypes

    def __assign_one_by_cluster(candidate_clusters: List[PocketCluster], datapoints: List[AtlasDatapoint]) -> None:
        """
        Uses the Hungarian algorithm to assign one datapoint to one cluster each. This is crucial, if one pocket
        includes several residues of the same type and even more of the same type in one single set.
        Example:
        For Ligand Tyr there is a pocket called GAA which includes Gly and two Ala. In one structure a Tyr is
        interacting with three Alanines. Of these three two must be assigned to the two pocket clusters. The third one
        should be the most outlying one and is discarded.

        :param candidate_clusters: A list of the Clusters the datapoints should be assigned to
        :param datapoints: A list of the datapoints
        :return: None
        """
        restypes = __get_restypes(datapoints)
        for restype, number in restypes.items():
            res_dp_list = [x for x in datapoints if x.binder_restype == restype]
            dist_matrix: List[List[float]] = [[] for _ in range(len(res_dp_list))]
            for i, res_dp in enumerate(res_dp_list):
                for c in candidate_clusters:
                    if c.binder_restype == restype:
                        dist_matrix[i].append(
                            c.get_dist_to_centroid(res_dp, distance_factor, orient_factor, secor_factor))
            if number != itemset.items[restype]:
                maximum = max(max(dist_matrix)) + 1.0
                for i, res_dp in enumerate(res_dp_list):
                    for dummy_clusters in range(number - itemset.items[restype]):
                        dist_matrix[i].append(maximum)
            assignment = linear_sum_assignment(dist_matrix)
            count = 0
            for c in candidate_clusters:
                if c.binder_restype == restype:
                    for i, a in enumerate(assignment[1]):
                        if a == count:
                            c.future_members.append(res_dp_list[assignment[0][i]])
                            break
                    count += 1


    clusters: List[PocketCluster] = []
    # create one cluster per item. Use representatives of a random pocket of the atlas as initial centroids.
    random_dp = choice(atlas.datapoints)
    random_pocket = [dp for dp in atlas.datapoints if dp.ligand_origin == random_dp.ligand_origin]
    for binder_restype, count in itemset.items.items():
        dps = sample(list([dp for dp in random_pocket if dp.binder_restype == binder_restype]), count)
        for dp in dps:
            clusters.append(PocketCluster(dp.ligand_restype, binder_restype, dp))
    # apply k-means clustering iteratively as long as the assignment of elements to clusters changes
    while any(c.has_changed for c in clusters):
        if one_per_cluster:
            # assign only "number-of-clusters" datapoints from each ligand residue to closest cluster
            for dps in atlas.group_by_ligand_origin().values():
                __assign_one_by_cluster(clusters, dps)
        else:
            # assign each datapoint to that cluster with the lowest distance to centroid
            for dp in atlas.datapoints:
                __assign_to_best_cluster(clusters, dp)
        # update the centroids and has_changed of every cluster
        for c in clusters:
            c.update(distance_factor, orient_factor, secor_factor)
    if one_per_cluster and more_than_one_ok:
        for dp in atlas.datapoints:
            __assign_to_best_cluster(clusters, dp)
        for c in clusters:
            c.update(distance_factor, orient_factor, secor_factor)
    # postprocessing: Until the variance of any of the clusters no longer exceeds the threshold, remove the most
    # distant member from the cluster.
    for c in clusters:
        while c.variance > variance_threshold:
            lrm = c.get_least_representative_member(distance_factor, orient_factor, secor_factor)
            c.future_members = [m for m in c.members if m != lrm]
            c.update(distance_factor, orient_factor, secor_factor)
    return Pocket(atlas.datapoints[0].ligand_restype, itemset, support, clusters)


def mine_pockets(atlas: Atlas, max_pockets_per_ligand_restype: int = 5, min_itemset_cardinality: int = 2,
                 confidence_threshold: float = 0.05, support_threshold: float = 0.05, variance_threshold: float = 5.0,
                 distance_factor: float = 1.0, orient_factor: float = 2.0, secor_factor: float = 2.0,
                 n_workers: int = -1, itemset_observer: Callable[[int, int], None] = None,
                 clustering_observer: Callable[[int, int], None] = None, cardinality_base: float = 2.0
                 ) -> Dict[str, List[Pocket]]:
    """Facade function that extracts pocket item datasets from an atlas, performs a frequent itemset mining analysis
    specifically for every residue type, and applies clustering in order to identify pockets. These pockets are then
    returned, grouped by ligand residue type
    :param atlas: the collection of datapoints taken as basis for the analysis
    :param max_pockets_per_ligand_restype: the result will contain no more than this number of pockets per ligand
    residue type
    :param min_itemset_cardinality: frequent itemsets must have at least this cardinality. I.e., the resulting pockets
    will consists of at least this number of residues
    :param confidence_threshold: association rules that fall below this threshold are not applied
    :param support_threshold: itemsets below this value are not considered as frequent
    :param variance_threshold: if apply_clustering, the residues contained by every cluster are allowed to have at
    most this variance; if this is exceeded, the cluster is reduced in an iterative post-processing step
    :param distance_factor: weight for the distance between c_alpha positions in distance function
    :param orient_factor: weight for the difference between primary orientation in distance function
    :param secor_factor: weight for the difference between secondary orientation in distance function
    :param n_workers: number of parallel processes to exploit
    :param itemset_observer: watches the progress of itemset mining
    :param clustering_observer: watches the progress of clustering
    :param cardinality_base: used for fitness calculation of each itemset in 'support * sqrt(base)^|itemset|'.
    Higher values result in pockets with more items.
    :return: a dictionary that maps every ligand residue type to a collection of a pockets.
    confidence, and the filtered atlas.
    """
    if n_workers == -1:
        n_workers = os.cpu_count()

    def __filter_and_clusterize(_ligrestype: str, _itemset: PocketItemset, _support: float) -> Pocket:
        """
        worker function that filters and clusterizes the itemset based on the atlas specified in the parent function.
        :param _ligrestype:
        :param _itemset:
        :param _support:
        :return: the pocket calculated by the worker
        """
        _pocket_atlas = filter_atlas_by_itemset(atlas, _ligrestype, _itemset, False)
        _pocket = clusterize_pocket_atlas(_pocket_atlas, _itemset, _support, variance_threshold,
                                          distance_factor, orient_factor, secor_factor)
        return _pocket

    datasets = extract_pocket_item_datasets(atlas)
    lig_restypes: List[str] = []
    frequent_itemsets: List[PocketItemset] = []
    supports: List[float] = []
    datasets_items = datasets.items()
    total = len(datasets_items)
    processed = 0
    for lig_restype, dataset in datasets_items:
        dataset_itemsets: List[PocketItemset] = []
        items_added = 0
        for (support, itemset) in dataset.find_frequent_itemsets(confidence_threshold, support_threshold,
                                                                 min_itemset_cardinality, cardinality_base):
            if not any(itemset in i for i in dataset_itemsets):  # avoid that subsets of existing itemsets are returned
                lig_restypes.append(lig_restype)
                dataset_itemsets.append(itemset)
                supports.append(support)
                items_added += 1
                if items_added >= max_pockets_per_ligand_restype:
                    break
        processed += 1
        if itemset_observer is not None:
            itemset_observer(processed, total)
        frequent_itemsets.extend(dataset_itemsets)
    pockets: List[Pocket] = []
    total = len(lig_restypes)
    processed = 0
    with ThreadPoolExecutor(n_workers) as executor:
        for pocket in executor.map(__filter_and_clusterize, lig_restypes, frequent_itemsets, supports):
            pockets.append(pocket)
            processed += 1
            if clustering_observer is not None:
                clustering_observer(processed, total)
    result: Dict[str, List[Pocket]] = {}
    for pocket in pockets:
        if pocket.restype not in result:
            result[pocket.restype] = []
        result[pocket.restype].append(pocket)
    return result


class PocketsEncoder(AtlasEncoder):
    """
    Enables the json Encoder to also encode an PocketCluster, Pocket and PocketItemset
    """

    def default(self, obj):
        if isinstance(obj, PocketCluster):
            return {'ligand_restype': obj.ligand_restype,
                    'binder_restype': obj.binder_restype,
                    'members': obj.members,
                    'centroid': obj.centroid,
                    'variance': obj.variance,
                    'future_members': obj.future_members,
                    'has_changed': obj.has_changed}
        if isinstance(obj, PocketItemset):
            return {'items': obj.items}
        if isinstance(obj, Pocket):
            return {'restype': obj.restype,
                    'itemset': obj.itemset,
                    'support': obj.support,
                    'clusters': obj.clusters}
        else:
            return AtlasEncoder.default(self, obj)
