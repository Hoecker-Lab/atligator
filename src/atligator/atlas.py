"""This module contains definitions for the data structure 'atlas', a collection of relative positions and orientations
of ligand and binder residues, seen from the subjective perspective of the binder. The module also contains methods for
the extraction of atlas data from PDB files.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-02-28
"""
import os
from collections import OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor
from typing import List, Iterator, Dict, Callable

# import h5py
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import Vector

from atligator.pdb_util import find_ligands, get_icoor_from_pdb_residue, is_amino_acid_residue, get_cbeta_position, \
    get_residues_within_radius, get_path, canonical_amino_acids


class AtlasDatapoint:
    """An atlas datapoint emerges from an interaction detected between two residues part of a ligand and binder chain,
    respectively."""

    def __init__(self, ligand_restype: str, binder_restype: str, ligand_origin: str or None, binder_origin: str or None,
                 ligand_atoms: Dict[str, Vector], binder_atoms: Dict[str, Vector]):
        """
        Creates a new instance.
        :param ligand_restype: ligand residue type
        :param binder_restype: binder residue type
        :param ligand_origin: path of the ligand residue where the data stems from
        :param binder_origin: path of the binder residue
        :param ligand_atoms: positions of all sidechain atoms of the ligand residue
        :param binder_atoms: positions of all sidechain atoms of the binder residue
        """
        self.ligand_restype = ligand_restype
        self.binder_restype = binder_restype
        self.ligand_origin = ligand_origin
        self.binder_origin = binder_origin
        self.ligand_atoms = ligand_atoms
        self.binder_atoms = binder_atoms

    def binder_calpha(self):
        """
        :return: the CA atom of the binder
        """
        return self.binder_atoms['CA']

    def binder_cbeta(self):
        """
        :return: the CB atom of the binder
        """
        return self.binder_atoms['CB']

    def binder_co(self):
        """
        :return: the C atom of the binder
        """
        return self.binder_atoms['C']

    def binder_orientation(self) -> Vector:
        """
        :return: the vector between CA and CB atoms of the binder
        """
        return self.binder_cbeta() - self.binder_calpha()

    def binder_secondary_orientation(self) -> Vector:
        """
        :return: the vector between CA and C atoms of the binder
        """
        return self.binder_co() - self.binder_calpha()

    def binder_to_pdb_residue(self, res_id: int) -> Residue:
        """Converts the information contained in the binder atoms of this datapoint into a PDB compliant residue
        :param res_id: the residue ID the result will have
        :return: the corresponding PDB residue
        """
        res = Residue(('', res_id, ' '), self.binder_restype, '')
        for atom_name, atom_position in self.binder_atoms.items():
            if self.binder_restype == 'GLY' and atom_name == 'CB':
                continue
            atom = Atom(atom_name, atom_position, 0.0, 1.0, ' ', atom_name, 0, atom_name[0])
            res.add(atom)
        return res

    def ligand_to_pdb_residue(self, res_id: int) -> Residue:
        """Converts the information contained in the ligand atoms of this datapoint into a PDB compliant residue
        :param res_id: the residue ID the result will have
        :return: the corresponding PDB residue
        """
        res = Residue(('', res_id, ' '), self.ligand_restype, '')
        for atom_name, atom_position in self.ligand_atoms.items():
            if self.ligand_restype == 'GLY' and atom_name == 'CB':
                continue
            atom = Atom(atom_name, atom_position, 0.0, 1.0, ' ', atom_name, 0, atom_name[0])
            res.add(atom)
        return res


class Atlas:
    """an atlas is a collection of atlas residues."""

    def __init__(self, datapoints: List[AtlasDatapoint], features: Dict = None):
        """creates a new atlas instance, which is represented explicitly by its datapoints.
        :param datapoints: a sequence of residues represented in the internal coordinate system of the respective
        context ligand.
        """
        self.datapoints = datapoints
        if features is None:
            self.features = {}
        else:
            self.features = features

    def group_by_ligand_origin_dep(self) -> OrderedDict:  # OrderedDict[str, List[AtlasResidue]]
        """Disjointly and completely partitions the data points contained in this atlas by their ligand origin.
        Was discontinued because left behind same-ligand-origin datapoints that are not in a row.
        :return: maps the ligand origin string to the subset of atlas residues it is connected with.
        """
        result: OrderedDict = OrderedDict({})  # OrderedDict[str, List[AtlasResidue]]
        curr_lig_origin = None
        curr_datapoints = []
        for datapoint in self.datapoints:
            if datapoint.ligand_origin != curr_lig_origin:
                if len(curr_datapoints) > 0:
                    result[curr_lig_origin] = curr_datapoints
                    curr_datapoints = []
                curr_lig_origin = datapoint.ligand_origin
            curr_datapoints.append(datapoint)
        if len(curr_datapoints) > 0:
            result[curr_lig_origin] = curr_datapoints
        return result

    def group_by_ligand_origin(self) -> Dict:  # Dict[str, List[AtlasResidue]]
        """Disjointly and completely partitions the data points contained in this atlas by their ligand origin.
        :return: maps the ligand origin string to the subset of atlas residues it is connected with.
        """
        result = defaultdict(list)
        for datapoint in self.datapoints:
            result[datapoint.ligand_origin].append(datapoint)
        return result

    def filter(self, ligand_restype: str = None, binder_restype: str = None) -> 'Atlas':
        """Filters this atlas by a given ligand and binder residue type.
        :param ligand_restype: if None, no filtering by ligand is applied.
        :param binder_restype: if None, no filtering by binder is applied.
        :return: the filtered atlas.
        """
        filtered_dps = self.datapoints[:]
        if ligand_restype is not None:
            filtered_dps = [dp for dp in filtered_dps if dp.ligand_restype == ligand_restype]
        if binder_restype is not None:
            filtered_dps = [dp for dp in filtered_dps if dp.binder_restype == binder_restype]
        return Atlas(filtered_dps)

    def get_stats(self, ligand_per_binder: bool = False) -> Dict[str, Dict[str, int]]:
        """Extracts statistical information about the distribution of datapoints from this atlas.
        :param ligand_per_binder: if True, the dictionary maps ligand to binder, else vice versa.
        :return: a dictionary that indicates the number of data points for every ligand (first key) and binder
        (second key) residue type.
        """
        stats: Dict[str, Dict[str, int]] = {}
        for dp in self.datapoints:
            key1 = dp.binder_restype if ligand_per_binder else dp.ligand_restype
            key2 = dp.ligand_restype if ligand_per_binder else dp.binder_restype
            if key1 not in stats:
                stats[key1] = {}
            if key2 not in stats[key1]:
                stats[key1][key2] = 0
            stats[key1][key2] = stats[key1][key2] + 1
        return stats

    def to_pdb_structure(self, structure_name: str) -> Structure:
        """
        :return: A PDB compliant structure that contains the datapoints of this atlas. By convention, chain A contains
        different conformations of the ligand residues of all datapoints, chains B contains binder residues.
        """
        model = Model(0)
        resid = 1
        ligand_chain = Chain('A')
        binder_chain = Chain('B')
        for dp in self.datapoints:
            ligres = dp.ligand_to_pdb_residue(resid)
            ligand_chain.add(ligres)
            resid += 1
        for dp in self.datapoints:
            binres = dp.binder_to_pdb_residue(resid)
            binder_chain.add(binres)
            resid += 1
        model.add(ligand_chain)
        model.add(binder_chain)
        struc = Structure(structure_name)
        struc.add(model)
        return struc

    def __len__(self):
        return len(self.datapoints)


def find_datapoints(model: Model, ir_default: float = 4.0, ir_hbond: float = 6.0, ir_aromatic: float = 6.0,
                    ir_ionic: float = 8.0, minlen: int = 4, maxlen: int = 20, include_hydrogens: bool = False,
                    alternative_lig_aa: str or None = None, restrict_to_alternative: bool = True,
                    allow_self_interactions: bool = False, skip_bb_atoms: bool = False) -> Iterator[AtlasDatapoint]:
    """Processes a PDB structure, identifies ligand chains in there, looks for interactions between ligand residues and
    binder residues. Converts each interacting binder residue into an atlas residue in internal coordinates, and yields
    all of these as results.
    :param model: the PDB model to process
    :param ir_default: interaction radius for interaction between arbitrary atoms in Angstrom
    :param ir_hbond: interaction radius for hydrogen bond interactions in Angstrom
    :param ir_aromatic: interaction radius for aromatic interactions in Angstrom
    :param ir_ionic: interaction radius for ionic interactions in Angstrom
    :param minlen: the minimum length of a chain in order to be considered as a ligand
    :param maxlen: the maximum length of a chain in order to be considered as a ligand
    :param include_hydrogens: whether hydrogen atoms shall be considered for interaction identification
    :param alternative_lig_aa: If given this alternative name of non-canonical aas will be also searched for as a
    ligand residue
    :param restrict_to_alternative: If True: If alternative Ligand amino acid is given restrict to this type.
    :param allow_self_interactions: If True interactions within one chain are taken into account (Slows down
    calculation heavily!!)
    :return: a generator of AtlasDatapoint instances
    """
    ligands = find_ligands(model, minlen, maxlen, alternative_lig_aa)
    for ligand in ligands:
        for ligand_residue in ligand.get_residues():
            # If alternative Ligand amino acid is given restrict to this type.
            if alternative_lig_aa is not None and \
                    restrict_to_alternative and not ligand_residue.get_resname() == alternative_lig_aa:
                continue
            if not (is_amino_acid_residue(ligand_residue, extended=alternative_lig_aa)):
                continue
            ligres: str = ligand_residue.get_resname()
            icoor = get_icoor_from_pdb_residue(ligand_residue)
            for binder_residue in get_residues_within_radius(model, ligand, ligand_residue, ir_default,
                                                             ir_hbond, ir_aromatic, ir_ionic,
                                                             include_hydrogens=include_hydrogens,
                                                             allow_ligand_residues=allow_self_interactions,
                                                             skip_bb_atoms=skip_bb_atoms):
                if not (is_amino_acid_residue(binder_residue)):
                    continue
                ligand_atoms: Dict[str, Vector] = {}
                for ligand_atom in ligand_residue.get_atoms():
                    if not ligand_atom.element == 'H':
                        ligand_atoms[ligand_atom.get_id()] = icoor.external_to_internal(ligand_atom.get_vector(), False)
                if 'CB' not in ligand_atoms:
                    ligand_atoms['CB'] = icoor.external_to_internal(get_cbeta_position(ligand_residue), False)
                binder_atoms: Dict[str, Vector] = {}
                for binder_atom in binder_residue.get_atoms():
                    if not binder_atom.element == 'H':
                        binder_atoms[binder_atom.get_id()] = icoor.external_to_internal(binder_atom.get_vector(), False)
                if 'CB' not in binder_atoms:
                    binder_atoms['CB'] = icoor.external_to_internal(get_cbeta_position(binder_residue), False)
                ligand_origin = get_path(ligand_residue)
                binder_origin = get_path(binder_residue)
                yield AtlasDatapoint(ligres, binder_residue.get_resname(), ligand_origin, binder_origin,
                                     ligand_atoms, binder_atoms)


def generate_atlas(filenames: List[str], ir_default: float = 4.0, ir_hbond: float = 6.0, ir_aromatic: float = 6.0,
                   ir_ionic: float = 8.0, minlen: int = 4, maxlen: int = 20, n_workers: int = -1,
                   include_hydrogens: bool = False, include_alternative_models: bool = False,
                   observer: Callable[[str], None] = None, alternative_lig_aa: str = None,
                   restrict_to_alternative: bool = True, allow_self_interactions: bool = False,
                   skip_bb_atoms: bool = True) -> Atlas:
    """Generates an atlas from a list of PDB files. In a multithreaded way, datapoints are collected and combined into
    a single atlas instance, which is returned as result.
    :param filenames: a list of PDB file names to process
    :param ir_default: interaction radius for interaction between arbitrary atoms in Angstrom
    :param ir_hbond: interaction radius for hydrogen bond interactions in Angstrom
    :param ir_aromatic: interaction radius for aromatic interactions in Angstrom
    :param ir_ionic: interaction radius for ionic interactions in Angstrom
    :param minlen: minimum chain length. see find_datapoints
    :param maxlen: maximum chain length. see find_datapoints
    :param n_workers: size of the threadpool, or -1 if it shall equal the number of physical CPU cores
    :param include_hydrogens: whether hydrogen atoms are to be considered. see find_datapoints
    :param include_alternative_models: whether to include data from models different from the first model of the
    structure (i.e., whether not to ignore alternative NMR conformations)
    :param observer: a callback function that receives the file name of every processed input structure
    :param alternative_lig_aa: If given this alternative name of non-canonical aas will be also searched for as a
    ligand residue
    :param restrict_to_alternative: If True: If alternative Ligand amino acid is given restrict to this type.
    :param allow_self_interactions: If True interactions within one chain are taken into account (Slows down
    calculation heavily!!)
    :return: The atlas created from the PDB files based on the given parameters
    """
    if n_workers == -1:
        n_workers = os.cpu_count()

    def __open_and_add_datapoints(_filename: str) -> List[AtlasDatapoint]:
        """Worker function that creates a PDB parser, opens a PDB file, calls find_datapoints, and returns the result.
        :param _filename: the PDB file to process
        :return: a list of datapoints that can be appended to an atlas
        """
        _parser = PDBParser(QUIET=True)
        # TODO pathlib
        _structure: Structure = _parser.get_structure(_filename.split('/')[-1], _filename)
        _datapoints: List[AtlasDatapoint] = []
        for _m in _structure:
            _datapoints.extend(find_datapoints(model=_m, ir_default=ir_default, ir_hbond=ir_hbond,
                                               ir_aromatic=ir_aromatic, ir_ionic=ir_ionic, minlen=minlen, maxlen=maxlen,
                                               include_hydrogens=include_hydrogens,
                                               alternative_lig_aa=alternative_lig_aa,
                                               restrict_to_alternative=restrict_to_alternative,
                                               allow_self_interactions=allow_self_interactions,
                                               skip_bb_atoms=skip_bb_atoms
                                               ))
            if not include_alternative_models:
                break
        if observer is not None:
            observer(_filename)
        return _datapoints

    atlas = Atlas([], features={"ir_default": ir_default, "ir_hbond": ir_hbond, "ir_aromatic": ir_aromatic,
                                "ir_ionic": ir_ionic, "minlen": minlen, "maxlen": maxlen,
                                "include_hydrogens": include_hydrogens,
                                "include_alternative_models": include_alternative_models,
                                "alternative_lig_aa": alternative_lig_aa,
                                "restrict_to_alternative": restrict_to_alternative,
                                "allow_self_interactions": allow_self_interactions, "skip_bb_atoms": skip_bb_atoms})
    with ThreadPoolExecutor(n_workers) as executor:
        for dps in executor.map(__open_and_add_datapoints, filenames):
            atlas.datapoints.extend(dps)
    return atlas
