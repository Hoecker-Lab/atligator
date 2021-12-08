"""This module contains data structures and operations for the handling of designed ligands and designed residues, based
on which the prediction of preferred residue types shall be undertaken.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-02-28
"""

from typing import List, Tuple, Dict

from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import Vector
from Bio.SeqUtils import seq1

from atligator.geometry import InternalCoordinates, Matrix
from atligator.pdb_util import is_amino_acid_residue, get_cbeta_position

# typedefs
LigandInfo = Tuple[str, str, str]  # (chain_id, residue_id, new_residue_type)
BinderInfo = Tuple[str, str]  # (chain_id, residue_id)
CouplingInfo = Tuple[str, List[str]]  # (chain_id, {residue_id})


class AbstractResidue:
    """A residue is defined by a dictionary of atoms in absolute coordinates."""

    def __init__(self, atoms: Dict[str, Vector]):
        """Creates a new residue.
        :param atoms
        """
        self.atoms = atoms

    def calpha(self) -> Vector:
        return self.atoms['CA']

    def cbeta(self) -> Vector:
        return self.atoms['CB']

    def co(self) -> Vector:
        return self.atoms['C']

    def orientation(self) -> Vector:
        """The orientation vector of a designed residue corresponds to the vector between C_alpha and C_beta atoms.
        :return: the orientation vector in absolute coordinates
        """
        return self.cbeta() - self.calpha()

    def secondary_orientation(self) -> Vector:
        """The secondary orientation vector of a designed residue corresponds to the vector between C_alpha and C_O
        atoms.
        :return: the secondary orientation vector in absolute coordinates
        """
        return self.co() - self.calpha()


def get_icoor_from_aresidue(residue: AbstractResidue) -> InternalCoordinates:
    """Defines an internal coordinate system relative to a given context residue as follows:
    - the origin is the context C_alpha atom
    - the x-axis is defined by the vector between context C_alpha and C_beta.
    - C_alpha, C_beta, and C lie within the xy-plane, such that the internal z axis denotes the normal of this plane.
    - C_N has always a negative z coordinate (when considering L amino acids)
    - lengths are preserved, no scaling is applied when compared to the external coordinate system.
    - C points to positive y direction! -> y is now -y, otherwise D-AAs are created
    :param residue: the context residue of the created internal coordinate system
    :return: a corresponding internal coordinate system
    """
    # basis vectors for the internal coordinate system
    internal_x = (residue.cbeta() - residue.calpha()).normalized()
    planar = (residue.co() - residue.calpha()).normalized()
    internal_z = (internal_x ** planar).normalized()
    internal_y = -(internal_x ** internal_z).normalized()
    # create the rotation matrix
    rotation_matrix: Matrix = [
        [internal_x[0], internal_x[1], internal_x[2]],
        [internal_y[0], internal_y[1], internal_y[2]],
        [internal_z[0], internal_z[1], internal_z[2]]]
    # return the internal coordinates. The translation vector equals the external C_alpha vector.
    return InternalCoordinates(residue.calpha(), rotation_matrix)


class DesignedLigandResidue(AbstractResidue):
    """A designed ligand residue is a residue with an individual position and orientation, which has in addition a
    ligand residue type assigned."""

    def __init__(self, atoms: Dict[str, Vector], chain: str, residue_id: int, restype: str):
        """Creates a designed ligand residue.
        :param atoms: the atoms of the residue and their positions
        :param chain: the identifier of the chain this residue is part of
        :param residue_id: a numerical residue identifier unique within the chain
        :param restype: the amino acid type of this designed ligand residues in three letter code
        """
        super().__init__(atoms)
        self.chain = chain
        self.residue_id = residue_id
        self.restype = restype


class DesignedLigand:
    """instances of this class represent designed ligands as a sequence of consecutive designed residues."""

    def __init__(self, residues: List[DesignedLigandResidue]):
        """Creates a new designed ligand instance.
        :param residues: designed residues comprising the initial contents of the ligand
        """
        self.residues = residues

    def get_title(self) -> str:
        """
        :return: a human readable representation of the ligand sequence.
        """
        title = ""
        for res in self.residues:
            title += seq1(res.restype)
        return title

    def get_min_dim(self) -> Tuple[int, int, int]:
        """
        :return: a tuple containing the minmal x, y, and z coordinates of any residues C_alpha atom
        """
        min_x: int = None
        min_y: int = None
        min_z: int = None
        for res in self.residues:
            if min_x is None or res.calpha()[0] < min_x:
                min_x = res.calpha()[0]
            if min_y is None or res.calpha()[1] < min_y:
                min_y = res.calpha()[1]
            if min_z is None or res.calpha()[2] < min_z:
                min_z = res.calpha()[2]
        return min_x, min_y, min_z

    def get_max_dim(self) -> Tuple[int, int, int]:
        """
        :return: a tuple containing the maximum x, y, and z coordinates of any residues C_alpha atom
        """
        max_x: int = None
        max_y: int = None
        max_z: int = None
        for res in self.residues:
            if max_x is None or res.calpha()[0] > max_x:
                max_x = res.calpha()[0]
            if max_y is None or res.calpha()[1] > max_y:
                max_y = res.calpha()[1]
            if max_z is None or res.calpha()[2] > max_z:
                max_z = res.calpha()[2]
        return max_x, max_y, max_z


def get_designed_ligand_from_scaffold(scaffold: Structure, ligand_info: List[LigandInfo]) -> DesignedLigand:
    """Extracts a designed ligand from a scaffold PDB structure and a list of ligand information structures.
    :param scaffold: a scaffold PDB structure that contains both the binder and the ligand in different chains
    :param ligand_info: describes how the designed ligand is extracted from the scaffold
    :return: a designed ligand constructed from the scaffold and the ligand info
    """
    residues: List[DesignedLigandResidue] = []
    for chain in scaffold.get_chains():
        for (chain_id, residue_id, new_residue_type) in ligand_info:
            if chain.get_id() == chain_id:
                for residue in chain.get_residues():
                    if str(residue.get_id()[1]) == residue_id:
                        if is_amino_acid_residue(residue):
                            atoms: Dict[str, Vector] = {}
                            for atom in residue.get_atoms():
                                atoms[atom.get_name()] = atom.get_vector()
                            if 'CB' not in atoms:
                                atoms['CB'] = get_cbeta_position(residue)
                            residues.append(DesignedLigandResidue(atoms, chain_id, int(residue_id), new_residue_type))
    return DesignedLigand(residues)


class DesignedBinderResidue(AbstractResidue):
    """A designed binder residue is a residue with an individual position and orientation, with no residue type
     information attached."""

    def __init__(self, atoms: Dict[str, Vector], chain: str, residue_id: int, original_residue_type: str):
        """Creates a new designed binder residue from existing data parameters.
        :param atoms: the atoms of the residue and their positions
        :param chain: the identifier of the chain this residue is part of
        :param residue_id: id of the residue in the scaffold's binder chain
        :param original_residue_type: residue type in the scaffold
        """
        super().__init__(atoms)
        self.chain = chain
        self.residue_id = residue_id
        self.original_residue_type = original_residue_type

    def transform_coor(self, icoor: InternalCoordinates, ext_to_int: bool) -> 'DesignedBinderResidue':
        """Converts this binder residue's atom coordinates from external to internal representation or vice versa.
        :param icoor: the coordinate system to use for transformation
        :param ext_to_int: if True, convert external to internal; if False, convert internal to external coordinates.
        :return: a transformed binder residue
        """
        trans_coor: Dict[str, Vector] = {}
        for name, pos in self.atoms.items():
            trans_pos = icoor.external_to_internal(pos, False) if ext_to_int else icoor.internal_to_external(pos, False)
            trans_coor[name] = trans_pos
        return DesignedBinderResidue(trans_coor, self.chain, self.residue_id, self.original_residue_type)


class DesignedBinder:
    """Instances of this class represent designed binders as a collection of (not necessarily connected) residues."""

    def __init__(self, residues: List[DesignedBinderResidue]):
        """Creates a new designed binder instance.
        :param residues initial collection of residues
        """
        self.residues = residues


def get_designed_binder_from_scaffold(scaffold: Structure, binder_info: List[BinderInfo]) -> DesignedBinder:
    """Extracts a designed binder from a scaffold PDB structure and a list of binder information structures.
    Adds a CB atom to the scaffold binder residues, if it is not there!!!

    :param scaffold: a scaffold PDB structure that contains both the binder and the ligand in different chains
    :param binder_info: explains how to extract designed binder residues from the scaffold
    :return: a designed binder with the extracted residues
    """
    residues: List[DesignedBinderResidue] = []
    for chain in scaffold.get_chains():
        for (chain_id, residue_id) in binder_info:
            if chain.get_id() == chain_id:
                for residue in chain.get_residues():
                    if str(residue.get_id()[1]) == residue_id:
                        original_restype = residue.get_resname()

                        # Add a dummy CB atom to the residues, if none is there. (Will show up in the structure!)
                        if 'CB' not in (atm.id for atm in residue.get_atoms()) and residue.resname != "GLY":
                            residue.add(Atom("CB", get_cbeta_position(residue), bfactor=0, occupancy=1.0,
                                             altloc=" ", fullname=f" CB ", serial_number=0, element="C"))
                        if is_amino_acid_residue(residue):
                            atoms: Dict[str, Vector] = {}
                            for atom in residue.get_atoms():
                                atoms[atom.get_name()] = atom.get_vector()
                            if residue.resname == "GLY":  # We need a dummy CB for Glycine
                                atoms["CB"] = get_cbeta_position(residue)
                            residues.append(DesignedBinderResidue(atoms, chain_id, int(residue_id), original_restype))
    return DesignedBinder(residues)
