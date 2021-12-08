"""This module contains utility functions for reading or manipulating entries of PDB files.

:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-05
"""

import gzip
import logging
import os
import urllib
from http.client import HTTPResponse
from math import degrees, sqrt
from shutil import copyfileobj
from typing import Iterator, List, Tuple, Dict, Union, TextIO

import Bio.PDB
import numpy as np
from Bio.PDB import Select, PDBList
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import Vector
from Bio.SeqUtils import seq1
from numpy import ndarray
from scipy import spatial

from atligator.geometry import InternalCoordinates, Matrix

# typedef "forwarded" from module atligator.structure
LigandInfo = Tuple[str, str, str]  # (chain_id, residue_id, new_residue_type)
BinderInfo = Tuple[str, str]  # (chain_id, residue_id)

canonical_amino_acids: List[str] = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
canonical_amino_acids_short: List[str] = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                                          'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
aa_3to1_conv: Dict[str, str] = dict(zip(canonical_amino_acids, canonical_amino_acids_short))
aa_1to3_conv: Dict[str, str] = dict(zip(canonical_amino_acids_short, canonical_amino_acids))

# TODO add alternative names that should be treated as aas:
# https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py


# source: http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
# Methionine: Biswal et al.: Sulfuric hydrogen bonds
# TODO pi ion interactions...
sidechain_functions: Dict[Tuple[str, str], List[str]] = {
    ('ARG', 'NE'): ['don'], ('ARG', 'NH1'): ['pos', 'don'], ('ARG', 'NH2'): ['pos', 'don'],
    ('ASN', 'ND2'): ['don'], ('ASN', 'OD1'): ['acc'],
    ('ASP', 'OD1'): ['neg', 'acc'], ('ASP', 'OD2'): ['neg', 'acc'],
    ('GLN', 'NE2'): ['don'], ('GLN', 'OE1'): ['acc'],
    ('GLU', 'OE1'): ['neg', 'acc'], ('GLU', 'OE2'): ['neg', 'acc'],
    ('HIS', 'CG'): ['pi'], ('HIS', 'CD2'): ['pi'], ('HIS', 'CE1'): ['pi'],
    ('HIS', 'ND1'): ['pos', 'don', 'acc', 'pi'], ('HIS', 'NE2'): ['pos', 'don', 'acc', 'pi'],
    ('LYS', 'NZ'): ['pos', 'don'],
    ('MET', 'SD'): ['acc'],
    ('PHE', 'CG'): ['pi'], ('PHE', 'CD1'): ['pi'], ('PHE', 'CD2'): ['pi'],
    ('PHE', 'CE1'): ['pi'], ('PHE', 'CE2'): ['pi'], ('PHE', 'CZ'): ['pi'],
    ('SER', 'OG'): ['don', 'acc'],
    ('THR', 'OG1'): ['don', 'acc'],
    ('TRP', 'NE1'): ['don', 'pi'], ('TRP', 'CG'): ['pi'], ('TRP', 'CD1'): ['pi'], ('TRP', 'CD2'): ['pi'],
    ('TRP', 'CE2'): ['pi'], ('TRP', 'CE3'): ['pi'], ('TRP', 'CZ2'): ['pi'], ('TRP', 'CZ3'): ['pi'],
    ('TRP', 'CH2'): ['pi'],
    ('TYR', 'OH'): ['don', 'acc'], ('TYR', 'CG'): ['pi'], ('TYR', 'CD1'): ['pi'], ('TYR', 'CD2'): ['pi'],
    ('TYR', 'CE1'): ['pi'], ('TYR', 'CE2'): ['pi'], ('TYR', 'CZ'): ['pi'],
    ('ALA', 'O'): ['acc'], ('ALA', 'N'): ['don'], ('ARG', 'O'): ['acc'], ('ARG', 'N'): ['don'],
    ('ASN', 'O'): ['acc'], ('ASN', 'N'): ['don'], ('ASP', 'O'): ['acc'], ('ASP', 'N'): ['don'],
    ('CYS', 'O'): ['acc'], ('CYS', 'N'): ['don'], ('GLN', 'O'): ['acc'], ('GLN', 'N'): ['don'],
    ('GLU', 'O'): ['acc'], ('GLU', 'N'): ['don'], ('GLY', 'O'): ['acc'], ('GLY', 'N'): ['don'],
    ('HIS', 'O'): ['acc'], ('HIS', 'N'): ['don'], ('ILE', 'O'): ['acc'], ('ILE', 'N'): ['don'],
    ('LEU', 'O'): ['acc'], ('LEU', 'N'): ['don'], ('LYS', 'O'): ['acc'], ('LYS', 'N'): ['don'],
    ('MET', 'O'): ['acc'], ('MET', 'N'): ['don'], ('PHE', 'O'): ['acc'], ('PHE', 'N'): ['don'],
    ('PRO', 'O'): ['acc'], ('PRO', 'N'): ['don'], ('SER', 'O'): ['acc'], ('SER', 'N'): ['don'],
    ('THR', 'O'): ['acc'], ('THR', 'N'): ['don'], ('TRP', 'O'): ['acc'], ('TRP', 'N'): ['don'],
    ('TYR', 'O'): ['acc'], ('TYR', 'N'): ['don'], ('VAL', 'O'): ['acc'], ('VAL', 'N'): ['don'],
}

ABC = ("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
       "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
       "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
       "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
       "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

# Conversion dictionaries for chain processing within ALARMS
chain_conv = dict(zip(ABC, range(len(ABC))))
chain_inv_conv = dict(zip(range(len(ABC)), ABC))


def get_interaction_types(restype1: str, atom1: str, restype2: str, atom2: str) -> List[str or None]:
    """Determines a specific interaction type between two atoms of different residues.
    :param restype1: three letter code of first residue
    :param atom1: atom code for first atom
    :param restype2: three letter code for second residue
    :param atom2: atom code for second atom
    :return: one of ['ionic', 'hbond', 'aromatic', 'unspecific']
    """
    result = []
    try:
        sf1 = sidechain_functions[(restype1, atom1)]
        sf2 = sidechain_functions[(restype2, atom2)]

        if 'pos' in sf1 and 'neg' in sf2:
            result.append('ionic')
        elif 'neg' in sf1 and 'pos' in sf2:
            result.append('ionic')

        if 'don' in sf1 and 'acc' in sf2:
            result.append('hbond')
        elif 'acc' in sf1 and 'don' in sf2:
            result.append('hbond')

        if 'pi' in sf1 and 'pi' in sf2:
            result.append('aromatic')
    except KeyError:
        pass
    result.append(None)

    return result


def _get_child_by_id(parent: Entity, child_id: str) -> Entity or None:
    """
    returns the child entity with specified ID in the PDB parent entity.
    This is not supposed to be executed by itself, rather wrapped by a function to check for the correct entity type.
    Thus, not tested yet.
    :param parent: the PDB entity to search.
    :param child_id: the ID of the child to find
    :return: the retrieved child, or None if not exists
    """
    try:
        return parent.child_dict[child_id]
    except KeyError:
        return None


def get_chain_by_id(model: Model, chain_id: str) -> Chain or None:
    """returns the chain with specified ID in the PDB structure.
    :param model: the PDB model to search
    :param chain_id: the ID of the chain to find
    :return: the retrieved chain, or None if not exists
    """
    if isinstance(model, Model):
        return _get_child_by_id(model, chain_id)
    else:
        raise TypeError


def get_path(residue: Residue) -> str:
    """
    Obtains a unique path name for a PDB residue. Example: '5XYZ/0/A/GLU250'
    :param residue:
    :return: a globally unique path residue
    """
    structure, model, chain, id_ = residue.get_full_id()
    return "/".join((structure, str(model), chain, residue.get_resname() + str(id_[1])))


def get_matching_sequences(model: Model, minlen: int, maxlen: int) -> Iterator[str]:
    """identifies chains matching the specified range in the structure and returns their sequences in IUPAC one letter
    code.
    :param model: the PDB model to search
    :param minlen: lower bound of the search range (inclusive)
    :param maxlen: upper bound of the search range (inclusive)
    :return: a generator for sequential strings
    """
    for chain in model.get_chains():
        s = ""
        aa_residues: List[Residue] = [r for r in chain.get_residues() if is_amino_acid_residue(r)]
        for res in aa_residues:
            s += seq1(res.get_resname())
        length = len(s)
        if minlen <= length <= maxlen:
            yield s


def is_canonical_amino_acid(resname: str) -> bool:
    """
    Determines whether the given name is part of the canonical amino acids
    :param resname: the amino acid type to look for
    :return: True if amino acid is accepted
    """
    return resname in canonical_amino_acids


def has_amino_acid_atoms(residue: Residue, glycine: bool = False) -> bool:
    """
    Determines whether the given residue includes C and CA and N atoms, as well as CB if it's not a glycine.
    :param residue: Bio.PDB residue entity to look for the atoms
    :param glycine: If true CB is not mandatory to pass
    :return: True if atoms are found
    """
    if glycine:
        return 'CA' in residue.child_dict and 'C' in residue.child_dict and 'N' in residue.child_dict
    else:
        return 'CA' in residue.child_dict and 'C' in residue.child_dict and 'N' in residue.child_dict and \
               'CB' in residue.child_dict


def is_amino_acid_residue(residue: Residue, extended: str or None = None) -> bool:
    """determines whether the given is an accepted amino acid residue, i.e.,
    it is part of the 20 canonical amino acids or an amino acid type specified in the extended parameter AND
    it has C and CA and N atoms, as well as CB if it's not a glycine.
    :param residue: the residue to check
    :param extended: if filled with str this string is an accepted aa name
    :return: whether the required properties are fulfilled for residue
    """
    resname = residue.get_resname()

    if is_canonical_amino_acid(resname):
        return has_amino_acid_atoms(residue, glycine=True if resname == "GLY" else False)
    elif extended is not None:
        return resname == extended and has_amino_acid_atoms(residue, glycine=True if resname == "GLY" else False)
    return False


def find_ligands(model: Model, minlen: int, maxlen: int, alternative_lig_aa: str or None = None) -> Iterator[Chain]:
    """returns a list of chains that satisfy the ligand criteria. As a first approximation, we assume that a ligand is
    a chain in the structure that contains >=minlen and <=maxlen residues.
    :param model: where to identify ligands
    :param minlen: minimum length of a chain for being a ligand
    :param maxlen: maximum length of a chain for being a ligand
    :param alternative_lig_aa: If given this alternative name of non-canonical aas will be also searched for as a
    ligand residue
    :return: a generator for chains being considered as ligands
    """
    for chain in model.get_chains():
        aa_residues: List[Residue] = [r for r in chain.get_residues()
                                      if is_amino_acid_residue(r, extended=alternative_lig_aa)]
        if minlen <= len(aa_residues) <= maxlen:
            yield chain


class NoBinderAtoms(Exception):
    """raised if no resdiue atoms are found in the binder chain."""

    def __init__(self, message="The binder chain does not contain suitable atoms!"):
        self.message = message
        super().__init__(self.message)


class NoLigandAtoms(Exception):
    """raised if no no atoms are found in other chains then the binder chain (ligand)."""

    def __init__(self, message="The chains other than the binder chain do not include suitable atoms!"):
        self.message = message
        super().__init__(self.message)


class NoResidueInDistance(Exception):
    """raised if no ligand residue is found in distance of the binder."""

    def __init__(self, message="No ligand residue is found in given distance of binder chain!"):
        self.message = message
        super().__init__(self.message)


def get_nonbinder_residues_within_radius(model: Model, binder_chain: Chain,
                                         ir_max: float = 8.0, ligand_min_len: int = 3,
                                         segmentation: bool = True,
                                         include_hydrogens: bool = True,
                                         check_canonical_aa: bool = False) -> Iterator[Residue]:
    """retrieves residues that are not part of the binder and lie within the interaction radius with any atom of a
    specified binder. Non-amino acid residues are ignored; optionally, hydrogen atoms are ignored.
    :param model: the PDB model from where to retrieve interactions
    :param binder_chain: residues of this chain are ignored as no self-interactions are recorded
    :param ir_max: interaction radius for interaction between arbitrary atoms in Angstrom
    :param segmentation: If True the binder residues get segmented. This allows to use min length of binder segments.
    :param ligand_min_len: Minimum length of segment of ligand chain
    :param include_hydrogens: whether hydrogen atoms shall be considered for all residues
    :param check_canonical_aa: whether to check if included residues are canonical amino acids or not
    :return: a generator of residues identified as interacting partners
    """

    # Define binder atoms: In binder chain and no hydrogens (default)
    binder_atoms = []
    if include_hydrogens:
        for at in binder_chain.get_atoms():
            if at.get_parent().get_id()[0] == " ":
                binder_atoms.append(at.get_coord())
    elif check_canonical_aa:  # checks if resname is canonical and has backbone atoms
        for at in binder_chain.get_atoms():
            if at.element != 'H' and is_amino_acid_residue(at.get_parent()):
                binder_atoms.append(at.get_coord())
    else:  # only checks if hetflag is " "
        for at in binder_chain.get_atoms():
            if at.element != 'H' and at.get_parent().get_id()[0] == " ":
                binder_atoms.append(at.get_coord())
    try:
        binder_atoms = np.concatenate(binder_atoms, axis=0).reshape((-1, 3))
    except ValueError:
        return NoBinderAtoms

    ligand_atoms = []
    ligand_atoms_res = []
    for ligand_ch in model.get_chains():
        chain_id = ligand_ch.get_id()
        chain_id_int = chain_conv[chain_id]
        if chain_id == binder_chain.get_id():
            continue
        for ligand_res in ligand_ch.get_residues():
            res_id = ligand_res.get_id()
            if res_id[0] != " " or (check_canonical_aa and is_amino_acid_residue(ligand_res)):
                continue
            for ligand_atom in ligand_res.get_atoms():
                if not include_hydrogens and ligand_atom.element == 'H':
                    continue

                # Select coordinates of corresponding atoms
                ligand_atoms.append(ligand_atom.get_coord())

                # Select residue id and chain id (converted to integer)
                ligand_atoms_res.append(res_id[1])
                ligand_atoms_res.append(chain_id_int)
    try:
        ligand_atoms = np.concatenate(ligand_atoms, axis=0).reshape((-1, 3))
    except ValueError:
        raise NoLigandAtoms
    ligand_atoms_res = np.asarray(ligand_atoms_res).reshape((-1, 2))

    # Build KD tree to find smallest distance to any ligand atom
    mytree = spatial.cKDTree(binder_atoms)
    dist, indexes = mytree.query(ligand_atoms)

    # Select those with distance below threshold and save on residue level
    res_in_dist = ligand_atoms_res[dist <= ir_max]
    if res_in_dist.size < 1:
        raise NoResidueInDistance
    # TODO warning if not possible?
    res_in_dist_uni = np.unique(res_in_dist, axis=0)  # Needs numpy >=1.13. If conda doesnt work (depend.) use pip.

    def consecutive(data, stepsize=1) -> List:
        """
        Returns a list of arrays. These arrays contain consecutive residues.
        :param data: Residues
        :param stepsize: to next residues. Should stay one.
        :return: List of Arrays with residues
        """
        return np.split(data, np.where(np.diff(data.flatten()[0::2].astype(int)) != stepsize)[0] + 1)

    # For segmentation
    if segmentation:
        for chain_int in chain_conv.values():
            chain_res_in_dist = np.take(res_in_dist_uni,
                                        np.where(res_in_dist_uni.flatten()[1::2] == chain_int)[0],
                                        axis=0)
            chain_res_in_dist = np.sort(chain_res_in_dist, axis=0)
            if chain_res_in_dist.size > 0:
                chain_res_in_dist = consecutive(chain_res_in_dist)
                for segment in chain_res_in_dist:
                    if segment.size >= ligand_min_len * 2:  # ligand minimal length!!!
                        for res in segment:
                            try:
                                yield model[chain_inv_conv[res[1]]][int(res[0])]
                            except KeyError:
                                pass
    # alternative mode without segmentation:
    else:
        for segment in res_in_dist_uni:
            yield model[chain_inv_conv[segment[1]]][int(segment[0])]


def get_residues_within_radius(model: Model, ligand: Chain, ligand_residue: Residue,
                               ir_default: float = 4.0, ir_hbond: float = 6.0,
                               ir_aromatic: float = 6.0, ir_ionic: float = 8.0,
                               include_hydrogens: bool = True,
                               allow_ligand_residues: bool = False, skip_bb_atoms: bool = False) -> Iterator[Residue]:
    """retrieves residues that are not part of the ligand (if allow_ligand_residues False) and lie within the
    interaction radius with any atom of a specified ligand. Non-amino acid residues are ignored; optionally
    for side-chain interactions, hydrogen atoms are ignored.
    :param model: the PDB model from where to retrieve interactions
    :param ligand: residues of this chain are ignored as no self-interactions are recorded
    :param ligand_residue: the ligand residue which is searched for interactions with non-ligand residues
    :param ir_default: interaction radius for interaction between arbitrary atoms in Angstrom
    :param ir_hbond: interaction radius for hydrogen bond interactions in Angstrom
    :param ir_aromatic: interaction radius for aromatic interactions in Angstrom
    :param ir_ionic: interaction radius for ionic interactions in Angstrom
    :param include_hydrogens: whether hydrogen atoms shall be considered for all residues
    :param allow_ligand_residues: If True also ligand_residues are allowed (not neighbours +- 3 residues!)
    :return: a generator of residues identified as interacting partners
    """
    if not is_amino_acid_residue(ligand_residue):
        pass  # TODO ?
    max_radius = max(ir_ionic, ir_aromatic, ir_hbond, ir_default)
    for chain in model.get_chains():

        # iteration over chains
        if chain.get_id() == ligand.get_id() and not allow_ligand_residues:
            continue
        forbidden = set()
        for binder_residue in chain:

            # iteration over residues in chain
            if allow_ligand_residues:
                if -3 < ligand_residue.get_id()[1] - binder_residue.get_id()[1] < 3:
                    continue

            added = False  # switch to see if residue is already found -> other atoms don't have to be iterated
            if binder_residue.get_id()[1] in forbidden:
                continue
            for binder_atom in binder_residue:

                # iteration over atoms in residue

                # skipping hydrogens or backbone atoms depending on arguments
                if not include_hydrogens and binder_atom.element == 'H':
                    continue
                if skip_bb_atoms and binder_atom.get_id() in ("N", "C", "O"):
                    continue

                for ligand_atom in ligand_residue:

                    # iteration over atoms in ligand residue

                    # skipping hydrogens or backbone atoms depending on arguments
                    if not include_hydrogens and ligand_atom.element == 'H':
                        continue
                    if skip_bb_atoms and ligand_atom.get_id() in ("N", "C", "O"):
                        continue

                    # requesting interaction type of specific atoms
                    itypes = get_interaction_types(ligand_residue.get_resname(), ligand_atom.get_id(),
                                                   binder_residue.get_resname(), binder_atom.get_id())
                    if len(itypes) == 1:
                        radius = ir_default
                    else:
                        radii = []
                        if 'ionic' in itypes:
                            radii.append(ir_ionic)
                        if 'aromatic' in itypes:
                            radii.append(ir_aromatic)
                        if 'hbond' in itypes:
                            radii.append(ir_hbond)
                        radius = max(radii)
                    distance: float = (binder_atom.get_vector() - ligand_atom.get_vector()).norm()
                    # Forbidden residues are neighbouring residues that will definitely be outside of maximum radius
                    # from the current ligand residue.
                    # 28 A are representing 8 A (as the the maximum distance from the most distant atoms of a residue)
                    # for the current ligand residue, the current binder_chain residue, the neighbour of current binder_chain
                    # residue as well as 4 A for the distance between two C-alpha atoms of neighbouring residues.
                    if distance > 28.0 + max_radius:
                        forbidden_neighbours = int((distance - 28.0) / 4)
                        curr_id = binder_residue.get_id()[1]
                        for i in range(forbidden_neighbours):
                            forbidden.add(curr_id + i)
                            forbidden.add(curr_id - i)
                    if distance < radius:
                        yield binder_residue
                        added = True
                        break
                if added:
                    break


def get_cbeta_position(residue: Residue) -> Vector:
    """returns the C_beta vector of this amino acid residue. For non-Gly residues, the real CB is returned. For Gly,
    this is simulated by geometric considerations.
    :param residue: the context residue
    :return: a vector (in absolute coordinates) representing the real or hypothetical C_beta location.
    """
    # 'real' C_beta exists
    if residue.has_id('CB'):
        return residue['CB'].get_vector()
    # C_beta does not exist; simulate it assuming that C, N, H_alpha and the virtual C_beta form a perfect tetrahedron
    # with C_alpha as centroid.
    calpha: Vector = residue['CA'].get_vector()
    co: Vector = residue['C'].get_vector()
    n_: Vector = residue['N'].get_vector()
    cent_co_n: Vector = (co + n_) ** 0.5
    cent_cb_h: Vector = cent_co_n + (calpha - cent_co_n) ** 2.0
    orth1 = (co - calpha).normalized()
    orth2 = (n_ - calpha).normalized()
    dir_cb_h = orth1 ** orth2
    length: float = (cent_co_n - co).norm()
    cbeta = cent_cb_h - dir_cb_h ** length
    return cbeta


def get_icoor_from_pdb_residue(residue: Residue) -> InternalCoordinates:
    """Defines an internal coordinate system relative to a given context PDB residue as follows:
    - the origin is the context C_alpha atom
    - the x-axis is defined by the vector between context C_alpha and C_beta.
    - C_alpha, C_beta, and C lie within the xy-plane, such that the internal z axis denotes the normal of this plane.
    - N has always a negative z coordinate (when considering L amino acids)
    - lengths are preserved, no scaling is applied when compared to the external coordinate system.
    - C points to positive y direction! -> y is now -y, otherwise D-AAs are created
    :param residue: the context residue for the internal coordinate system
    :return: an internal coordinate system with the described properties
    """
    calpha: Vector = residue['CA'].get_vector()
    cbeta: Vector = get_cbeta_position(residue)
    co: Vector = residue['C'].get_vector()
    # basis vectors for the internal coordinate system
    internal_x = (cbeta - calpha).normalized()
    planar = (co - calpha).normalized()
    internal_z = (internal_x ** planar).normalized()
    internal_y = -(internal_x ** internal_z).normalized()
    # create the rotation matrix
    rotation_matrix: Matrix = [
        [internal_x[0], internal_x[1], internal_x[2]],
        [internal_y[0], internal_y[1], internal_y[2]],
        [internal_z[0], internal_z[1], internal_z[2]]]
    # return the internal coordinates. The translation vector equals the external C_alpha vector.
    return InternalCoordinates(calpha, rotation_matrix)


def mutate_residue(residue: Residue, restype: str) -> None:
    """Applies a point mutation to a given residue. All atoms except for 'C', 'CA', 'N', 'O' are removed and the
    residue type is set according to the given value.
    :param residue: the residue to mutate
    :param restype: the new residue type
    """
    for scatom in [atom for atom in residue if atom.id not in ['C', 'CA', 'N', 'O']]:
        residue.detach_child(scatom.id)
    residue.resname = restype


def cut_termini(structure: Structure, chain_id: str = "A", chain_id_lig: str = "D",
                cut_off_binder: Tuple[int, int] = (52, 250), cut_off_ligand: Tuple[int, int] = (2, 9),
                rep_size: int = 42, delete_only: bool = False) -> None:
    """
    Trying to shorten an armadillo structure by cutting terminal repeat units. 5AEI = 42 res per unit.
    -> 7 units - 2 = 5 units. Starting with 11 - 291 shortened to 53 - 249. Ligand shortend from 1-10 to 3-8.
    To consider: terminal parts don't share same sequence.
    -> To go: copy aa sequence from 11-52 to 53-94 and from 250-291 to 208-249
              delete 11-52 and 250-291 and lig 1, 2, 9, 10.
    :param structure:
    :param chain_id:
    :param chain_id_lig:
    :param cut_off_binder:
    :param cut_off_ligand:
    :param rep_size:
    :param delete_only:
    :return: None
    """
    lower_end = cut_off_binder[0]
    upper_end = cut_off_binder[1]
    assert lower_end <= upper_end
    if not delete_only:
        for residue in structure[0][chain_id]:
            res_id: int = residue.get_id()[1]
            if res_id <= lower_end:
                mutate_residue(structure[0][chain_id][res_id + rep_size], residue.get_resname())
            elif res_id >= upper_end:
                mutate_residue(structure[0][chain_id][res_id - rep_size], residue.get_resname())
    to_remove = []
    for residue in structure[0][chain_id]:
        res_id: int = residue.get_id()[1]
        if res_id <= lower_end or res_id >= upper_end:
            to_remove.append(residue)
    for i in to_remove:
        structure[0][chain_id].detach_child(i.get_id())
    lower_end = cut_off_ligand[0]
    upper_end = cut_off_ligand[1]
    assert lower_end <= upper_end
    to_remove = []
    for residue in structure[0][chain_id_lig]:
        res_id: int = residue.get_id()[1]
        if res_id <= lower_end or res_id >= upper_end:
            to_remove.append(residue)
    for i in to_remove:
        structure[0][chain_id_lig].detach_child(i.get_id())


def replace_restypes(structure_file: str, old_type: str, new_type: str) -> None:
    """Replaces old_type residues with new_type residues everywhere in the provided PDB structure
    :param structure_file: PDB filename
    :param old_type: all occurrences of this residue type are to be replaced by new_type
    :param new_type
    """
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(structure_file, structure_file)
    for residue in structure.get_residues():
        if residue.resname == old_type:
            residue.resname = new_type
    io = PDBIO()
    io.set_structure(structure)
    io.save(structure_file)


class SecondaryScores:
    # TODO add rmsd to original from pymol_util
    def __init__(self, linearity: float, lin_p_res: float, phi: float, psi: float, phi_dev: float, psi_dev: float,
                 rmsd: float):
        self.linearity = linearity
        self.linearity_p_res = lin_p_res
        self.phi = phi
        self.psi = psi
        self.phi_dev = phi_dev
        self.psi_dev = psi_dev
        self.rmsd = rmsd


def n(vector: Tuple) -> ndarray:
    """
    Returning an numpy vector from a Tuple.
    :param vector: Tuple to convert to a vector
    :return: numpy vector
    """
    return np.array(vector)


def get_angle(a: Tuple, b: Tuple, c: Tuple) -> float or None:
    """
    Calculates an angle between the points a<->b<->c with b being the middle.
    Was tested with threedimensional points.
    :param a: point a with variables in a tuple
    :param b: point b with variables in a tuple
    :param c: point c with variables in a tuple
    :return: The angle in radians
    """
    a: ndarray = n(a)
    b: ndarray = n(b)
    c: ndarray = n(c)
    ba = a - b
    bc = c - b
    if np.linalg.norm(ba) == 0.0 or np.linalg.norm(bc) == 0.0:
        return None
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    cosine_angle = max((-1, min((1, cosine_angle))))
    angle: float = np.arccos(cosine_angle)
    return angle


def calc_dist_to_line(line_point_a: Tuple, line_point_b: Tuple, point_c: Tuple) -> float or None:
    """
    Calculates the distance of point c to a line crossing point a and b.
    :param line_point_a: coordinates of point a in a tuple.
    :param line_point_b: coordinates of point b in a tuple.
    :param point_c: coordinates of point c whose distance should be measured.
    :return: distance as a float
    """
    alpha = get_angle(line_point_a, line_point_b, point_c)
    if alpha is None:
        return None
    distance: float = np.sin(alpha) * np.linalg.norm((n(line_point_b) - n(point_c)))
    return distance


def get_binder_positions(structure: Structure, binder_info: List[BinderInfo]) -> List[int]:
    """
    Returns a list of absolute positions (starting at 1 for the first residue of the first chain) that match any
    content of the binder_info list (= residues of ligand binding site on binder_chain).
    Example: host consists of chain A, residue ids are 11, 12, 16, 17 of which 12 and 16 are part of the binding site.
    Ligand consists of residues 1, 2, 3 of chain D.
    The function would return [2, 3].
    :param structure: a Bio.PDB structure to identify binder_chain positions in
    :param binder_info: a list of tuples containing chain identifer and chain-relative residue id.
    :return: a list of absolute positions indicating the residues contained in binder_info
    """
    result: List[int] = []
    for i, res in enumerate(structure.get_residues()):
        for bi in binder_info:
            res: Residue
            if res.get_parent().get_id() == bi[0] and res.get_id()[1] == int(bi[1]):
                result.append(i + 1)
                break
    return result


def get_ligand_positions(structure: Structure, ligand_info: List[LigandInfo]) -> List[int]:
    """
    Returns a list of absolute positions (starting at 1 for the first residue of the first chain) that match any
    content of the ligand_info list.
    Example: host consists of chain A, residue ids are 11, 12, 16, 17. Ligand consists of residues 1, 2, 3 of chain D.
    The function would return [5, 6, 7].
    :param structure: a PDB structure to identify ligand positions in
    :param ligand_info: a list of tuples containing chain identifer and chain-relative residue id.
    :return: a list of absolute positions indicating the residues contained in ligand_info
    """
    result: List[int] = []
    for i, res in enumerate(structure.get_residues()):
        for li in ligand_info:
            res: Residue
            if res.get_parent().get_id() == li[0] and res.get_id()[1] == int(li[1]):
                result.append(i + 1)
                break
    return result


def inspect_pdb(pdb_id: str, ligand_info: List[LigandInfo], context_path: str = None) -> SecondaryScores:
    """
    inspects pdb file with name pdb_id.pdb in current directory and returns secondary
    scores like the linearity of the chain(s) other than A (-> peptide ligand) or the phi
    and psi angle averages of these chains.
    :param pdb_id: name of the pdb file without the '.pdb'
    :param ligand_info: contains descriptions of the residues of the ligand
    :param context_path: where it all happens
    :return: a SecondaryScore instance containing several attributes like phi/psi angle averages
    """
    for model in Bio.PDB.PDBParser().get_structure(pdb_id, pdb_id + '.pdb'):
        for chain in model:
            if chain.get_id() != ligand_info[0][0]:
                continue
            line_a: Tuple = chain[int(ligand_info[0][1])]['CA'].coord
            line_b: Tuple = chain[int(ligand_info[-1][1])]['CA'].coord
            linearity = 0.0
            n_res = len(chain)
            for residue in chain:
                for atom in residue:
                    if atom.get_id() in ('CA', 'C', 'N'):
                        if int(residue.get_id()[1]) not in (1, len(ligand_info)):
                            dist = calc_dist_to_line(line_a, line_b, atom.coord)
                            if dist is None:
                                continue
                            linearity += dist
            polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
            for poly in polypeptides:
                phipsi: List[List[float, float]] = poly.get_phi_psi_list()
                phps_deg: List[List[float]] = []
                for res in phipsi[1:-1]:
                    res_phps: List[float] = []
                    for i in range(2):
                        if res[i] is not None:
                            res_phps.append(degrees(res[i]))
                        else:
                            res_phps.append(0.0)
                    phps_deg.append(res_phps)
                phi_avg = 0.0
                psi_avg = 0.0
                tot = 0
                for i in range(len(phps_deg[1:-1])):
                    phi_avg += phps_deg[i][0]
                    psi_avg += phps_deg[i][1]
                    tot = i + 1
                phi_avg /= tot
                psi_avg /= tot
                phi_dev = 0
                psi_dev = 0
                for i in range(len(phps_deg[1:-1])):
                    phi_dev += (phps_deg[i][0] - phi_avg) ** 2
                    psi_dev += (phps_deg[i][1] - psi_avg) ** 2
                phi_dev = sqrt(phi_dev / tot)
                psi_dev = sqrt(psi_dev / tot)
                rmsd = get_rmsd(('' if context_path is None else context_path + os.sep) + 'scaffold',
                                pdb_id, len(ligand_info))
                return SecondaryScores(linearity, linearity / n_res, phi_avg, psi_avg, phi_dev, psi_dev, rmsd)


def get_rmsd(a_pdb: str, b_pdb: str, last_x: int = None) -> float:
    """
    Returns the rmsd value of the two pdb structures in parameters. If last_x parameter is given,
    Only the rmsd for the last 'last_x' residues is calculated. These last are also not taken into
    account for kabsch rotation-translation (very different positions would affect translation and rotation).

    :param a_pdb: Pdb id of first structure in curr. dir. '.pdb' will be added automatically.
    :param b_pdb: Pdb id of second structure in curr. dir. '.pdb' will be added automatically.
    :param last_x: Number of the last residues (peptide ligand?) to consider for RMSD
    :return: rmsd value
    """
    a = np.array(get_pdb_coors(a_pdb))
    b = np.array(get_pdb_coors(b_pdb))
    # Check if both have the same number of res
    assert len(a) == len(b)
    if last_x is None:
        last_x = len(b)
    b -= get_centroid(b[:-last_x])
    a -= get_centroid(a[:-last_x])
    a = rotate_kabsch(a, b, last_x)
    return calc_rmsd(a[-last_x:], b[-last_x:])


def get_pdb_coors(pdb_id: str) -> List[ndarray]:
    for model in Bio.PDB.PDBParser().get_structure(pdb_id, pdb_id + '.pdb'):
        coor_lst = []
        # TODO does it always start with chain A??
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_id() == 'CA':  # OR in ('CA', 'C', 'N'):
                        coor_lst.append(atom.coord)
        return coor_lst


def rotate_kabsch(a: ndarray, b: ndarray, not_include: int):
    r = get_u_kabsch(a[:-not_include], b[:-not_include])
    a = np.matmul(a, r)
    return a


def get_u_kabsch(a: ndarray, b: ndarray):
    """
    Returns the rotational matrix for rotation of a to fit b best. Calculated using
    kabsch algorithm.
    :param a: Matrix a including coordinates of first structure
    :param b: Matrix b including coordinates of second structure
    :return: Rotation matrix
    """
    c = np.matmul(np.transpose(a), b)
    # TODO check svd names and transpositons
    v, s, w = np.linalg.svd(c)
    d = (np.linalg.det(v) * np.linalg.det(w)) < 0.0

    if d:
        s[-1] = -s[-1]
        v[:, -1] = -v[:, -1]

    r = np.matmul(v, w)
    return r


def get_centroid(matrix: ndarray) -> ndarray:
    """
    Returning the x centroid values for a matrix with points in a x-dim. coordinate system.
    :param matrix: Matrix with several points and their coordinates
    :return: Centroid values for each dimension as a ndarray
    """
    centroid: ndarray = matrix.mean(axis=0)
    return centroid


def calc_rmsd(a: ndarray, b: ndarray) -> float:
    """
    Calculates the root mean square deviation for pairs of points with x coordinates (3 in 3D coor system) each.
    These points are given by ndarray matrices (a and b).
    :param a: Matrix a including coordinates of first structure
    :param b: Matrix b including coordinates of second structure
    :return: rmsd value
    """
    dimensions = len(a[0])
    n_points = len(a)
    rmsd = 0.0
    for a_i, b_i in zip(a, b):
        rmsd += sum([(a_i[i] - b_i[i]) ** 2.0 for i in range(dimensions)])
    return np.sqrt(rmsd / n_points)


def find_pdb_keys_in_scop_dir_cla(query: str, scop_dir_cla_file: Union[TextIO, HTTPResponse],
                                  decode: bool = False) -> Union[List[str], None]:
    """
    Parses a SCOP(e) directory classification file (dir.cla.scop*.txt) in order to find a list of PDB keys that match
    the given scop query (e.g., 'a', 'a.118.1', etc.)
    :param query: the class, fold, superfamily, or family to search
    :param scop_dir_cla_file: part of the SCOP(e) database that contains the classification
    :param decode: If True it tries to decode the lines in scop_dir_cla_file
    :return: list of PDB files that match the query, or None if the query was not well-formed
    """
    pdb_keys = set()
    query_split = query.split('.')
    if len(query_split) < 1:
        return None
    query_split = query_split if len(query_split[-1]) else query_split[:-1]
    for line in scop_dir_cla_file:
        if decode:
            line = line.decode("utf-8")
        entries = line.split()
        if not line.startswith('#') and len(entries) > 5:
            pdb_key = entries[1]
            scop_key = entries[3].split('.')
            add = True
            for i, sk in enumerate(scop_key):
                if i < len(query_split) and sk != query_split[i]:
                    add = False
                    break
            if add and pdb_key not in pdb_keys:
                pdb_keys.add(pdb_key)
    return sorted(list(pdb_keys))


def get_available_scop_dir_cla_file(scop_dir_cla_file: str = None) -> Tuple[Union[TextIO, HTTPResponse], bool]:
    """
    Returns the opened scop_dir_cla_file necessary for matching scop ids and pdb ids in a tuple with a boolean
    that shows if the file needs to be decoded or not.

    :param scop_dir_cla_file: file of the SCOP(e) database that contains the classification. Either string with path to
    the file locally or None. If None the file is downloaded from
    https://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.08-stable.txt
    :return: Tuple of the opened scop_dir_cla_file and a boolean that shows, if file needs to be decoded before use.
    """
    if scop_dir_cla_file is not None and os.path.isfile(scop_dir_cla_file):
        return open(scop_dir_cla_file, "r"), False
    logging.info("Downloading https://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.08-stable.txt to search for "
                 "matching patterns...")
    return urllib.request.urlopen("https://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.08-stable.txt"), True


def find_scop_i_from_pdbs(query: List[str], scop_dir_cla_file: str = None) \
        -> Dict[str, str] or None:
    """
    Parses a SCOP(e) directory classification file (dir.cla.scop*.txt) in order to find the scop identifiers for
    a list of PDB keys.
    :param query: a list of pdb ids you want the scop identifiers from
    :param scop_dir_cla_file: file of the SCOP(e) database that contains the classification. Either string with path to
    the file locally or None. If None the file is downloaded from
    https://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.08-stable.txt
    :return: list of PDB files that match the query, or None if the query was not well-formed
    """
    pdb_keys = set(query)
    pdb_scop_dict = dict((a, None) for a in query)
    with get_available_scop_dir_cla_file(scop_dir_cla_file=scop_dir_cla_file) as (fo, decode):
        for line in fo:
            if decode:
                line = line.decode("utf-8")
            entries = line.split()
            if not line.startswith('#') and len(entries) > 5:
                pdb_key = entries[1]
                if pdb_key in pdb_keys:
                    if pdb_scop_dict[pdb_key] is None:
                        scop_key = entries[3]
                        pdb_scop_dict[pdb_key] = [scop_key]
                    elif entries[3] not in pdb_scop_dict[pdb_key]:
                        scop_key = entries[3]
                        logging.info(f"Pdb file found twice: {pdb_key} in {pdb_scop_dict[pdb_key]} and {entries[3]}!")
                        pdb_scop_dict[pdb_key] += [scop_key]
    return pdb_scop_dict


def get_pdb_structure(pdb_key: str, dest_dir: str, local: bool = False, pdb_db_path: str = None) -> str or None:
    """
    Downloads a PDB file matching a given PDB key from a local or remote PDB database. Copies the file into a given
    destination directory and returns the path of the output file.
    If the pdb_key cannot be found in the local database (if local==True) searches in the global PDB database
    (via PDB.Bio.PDBList). If the corresponding pdb_key cannot be found in the database the returned path is not
    touched and can be checked if exists.

    :param pdb_key: The PDB key to find, e.g., 5AEI
    :param dest_dir: file will be copied into this folder as [pdb_key].pdb. If the path does not exist it will be
    created automatically.
    :param local: If True try to copy the file from local copy of pdb db (must be given in pdb_db_path and
    should rescolve : f"{pdb_db_path}/data/{pdb_key[1:3]}/pdb{pdb_key}.ent.gz"
    :param pdb_db_path: Directory that contains a copy of the PDB database. Must contain a folder 'data'.
    :return: the file name, if copying was successful
    """

    def return_remote():
        return PDBList().retrieve_pdb_file(pdb_key,
                                           pdir=dest_dir,
                                           file_format="pdb")

    if not local or pdb_db_path is None:
        pdb_path = return_remote()

    else:
        arch = f"{pdb_db_path}/data/{pdb_key[1:3]}/pdb{pdb_key}.ent.gz"
        pdb_path = f"{dest_dir}/{pdb_key}.pdb"
        if not os.path.exists(arch):
            pdb_path = return_remote()
        else:
            with gzip.open(arch, 'rb') as f_in, open(pdb_path, 'wb') as f_out:
                copyfileobj(f_in, f_out)
    return pdb_path


def get_structures_matching_scop_query(query: str, scop_dir_cla_file: str = None):
    """
    Finds all PDB files that match a given SCOP(e) query using a local copy of the SCOP(e) database
    (if scop_dir_cla_file is available) or

    :param query: the class, fold, superfamily, or family to search, e.g., 'a.118.1'
    :param scop_dir_cla_file: file of the SCOP(e) database that contains the classification. Either string with path to
    the file locally or None. If None the file is downloaded from
    https://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.08-stable.txt
    """
    if scop_dir_cla_file is not None and os.path.isfile(scop_dir_cla_file):
        with get_available_scop_dir_cla_file(scop_dir_cla_file=scop_dir_cla_file) as (scop_dir_cla_opened, decode):
            return find_pdb_keys_in_scop_dir_cla(query=query, scop_dir_cla_file=scop_dir_cla_opened, decode=decode)


def download_pdbs_by_scop_query(query: str, dest_dir: str, scop_dir_cla_file: str = None, local_db: bool = False,
                                pdb_db_path: str = None) -> None:
    """
    Finds all PDB files that match a given SCOP(e) query using a local copy of the SCOP(e) database, downloads the
    structures from a local copy of the PDB database, and stores the files in PDB format under the given destination
    directory.
    :param query: the class, fold, superfamily, or family to search, e.g., 'a.118.1'
    :param dest_dir: PDB structure files will be copied into this folder
    :param scop_dir_cla_file: part of the SCOP(e) database that contains the classification
    :param local_db: if True it searches for the file in a local copy of the protein database. Must be given in
    pdb_db_path
    :param pdb_db_path: Directory that contains a copy of the PDB database. Must contain a folder 'data'.
    """
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    pdbkeys = get_structures_matching_scop_query(query=query, scop_dir_cla_file=scop_dir_cla_file)
    print("Found no matching files." if not len(pdbkeys) else f"Found {len(pdbkeys)} matching pdb files.")
    for k in pdbkeys:
        out_file = get_pdb_structure(pdb_key=k, dest_dir=dest_dir, pdb_db_path=pdb_db_path, local=local_db)
        if not os.path.isfile(out_file):
            logging.warning(f"The file {out_file} cannot be found in the pdb database.")


class ChainSelect(Select):
    """
    A Select class to enable selecting the residues that should be written in an output structure file within
    Bio.PDB.
    """

    # TODO insert a 'remove waters and other non-aa-res ligands' option
    def __init__(self, allowed, chain_id="x"):
        super().__init__()
        self.chain_id = chain_id
        self.allowed = allowed

    def accept_residue(self, residue):
        """
        Only accept residues that are in self.allowed or that belong to the chain in defined self.chain_id (while being
        an amino acid residue)
        :param residue: the residue that get's tested
        :return: bool if this residue should be accepted or not
        """

        return residue in self.allowed or \
               (residue.get_parent().get_id() == self.chain_id and is_amino_acid_residue(residue))
