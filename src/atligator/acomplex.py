"""A ligand-binder complex consists of a designed ligand, a designed binder, and a sequence of aligned residues which
originate from an atlas. This module contains data structures and operations for creating and manipulating such
ligand-binder complexes.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-03-05
"""

from typing import List, Dict, Any

from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import Vector

from atligator.atlas import Atlas
from atligator.structure import AbstractResidue, DesignedLigand, DesignedBinder, LigandInfo, BinderInfo, CouplingInfo, \
    get_designed_binder_from_scaffold, get_designed_ligand_from_scaffold, get_icoor_from_aresidue


class ComplexDescriptor:
    """A complex descriptor complements a scaffold PDB file with additional information necessary for the construction
    of a designed ligand and binder out of the scaffold."""

    def __init__(self, ligand_info: List[LigandInfo], binder_info: List[BinderInfo], coupling_info: List[CouplingInfo]):
        """Creates a new complex descriptor instance from the data it is going to contain.
        :param ligand_info: a list of tuples, each containing chain id, residue id, and designed residue type
        :param binder_info: a list of tuples, each containing chain id and residue id
        :param coupling_info: a list of coupling information. The first element denotes the chain id, the second element
        is a list of coupled residues.
        """
        self.ligand_info = ligand_info
        self.binder_info = binder_info
        self.coupling_info = coupling_info


def read_complex_descriptor(filename: Any) -> ComplexDescriptor:
    """Reads in a descriptor file and returns a corresponding instance of ComplexDescriptor.
    :param filename: Pathname of a file to be opened and read or File itself
    :return a corresponding instance of ComplexDescriptor
    """

    def _get_properties(file):
        linfo: List[LigandInfo] = []
        binfo: List[BinderInfo] = []
        cinfo: List[CouplingInfo] = []
        stage = 0
        for line in file:
            if len(line.split()) <= 1:
                stage += 1
                continue
            content = line.split()
            if stage == 0:
                linfo.append((content[0], content[1], content[2]))
            elif stage == 1:
                binfo.append((content[0], content[1]))
            elif stage == 2:
                cinfo.append((content[0], content[1:]))
        return linfo, binfo, cinfo

    if not isinstance(filename, str):
        ligand_info, binder_info, coupling_info = _get_properties(filename)
    else:
        with open(filename, 'r') as fo:
            ligand_info, binder_info, coupling_info = _get_properties(fo)
    return ComplexDescriptor(ligand_info, binder_info, coupling_info)


class AlignedResidue(AbstractResidue):
    """An aligned residue is a residue with an individual position and orientation, which has in addition a binder
    residue type as well as a position of the originating residue in the designed ligand assigned."""

    def __init__(self, atoms: Dict[str, Vector], binder_restype: str, ligand_res_position: int, size_factor: int):
        """Creates a new instance.
        :param atoms: atoms part of this residues and their positions
        :param binder_restype: binder residue type inferred from the atlas
        :param ligand_res_position: the residue id relative to the ligand chain
        :param size_factor: number of datapoints of the original atlas
        """
        super().__init__(atoms)
        self.binder_restype = binder_restype
        self.ligand_res_position = ligand_res_position
        self.size_factor = size_factor


class LigandBinderComplex:
    """A ligand-binder complex is created from a scaffold structure. It consists of a designed ligand structure, a set
    of designed binder residues, and a set of aligned residues inferred from the atlas."""

    def __init__(self, ligand: DesignedLigand, binder: DesignedBinder, aligned_residues: List[AlignedResidue]):
        """Creates a new instance. All vectors are supposed to be in absolute coordinates (defined by the scaffold)
        :param ligand: the designed ligand part of the complex
        :param binder: the designed binder part of the complex
        :param aligned_residues: a set of aligned residue
        """
        self.ligand = ligand
        self.binder = binder
        self.aligned_residues = aligned_residues


def generate_complex(atlas: Atlas or None, scaffold: Structure, descriptor: ComplexDescriptor) -> LigandBinderComplex:
    """Generates a ligand-binder complex out of an atlas database, a scaffold structure, and a descriptor that refers
    to the scaffold structure.
    :param atlas: the atlas from which aligned residues will be inferred from
    :param scaffold: a scaffold PDB structure from which the position and orientation of residues of both the designed
    binder and the designed ligand will be extracted.
    :param descriptor: a descriptor file that contains information about ligand and binder residues and types
    :return: the generated ligand-binder complex
    """
    ligand = get_designed_ligand_from_scaffold(scaffold, descriptor.ligand_info)
    binder = get_designed_binder_from_scaffold(scaffold, descriptor.binder_info)
    aligned_residues: List[AlignedResidue] = []
    if atlas is not None:
        for lresid, ligres in enumerate(ligand.residues):
            icoor = get_icoor_from_aresidue(ligres)
            relevant_datapoints = [d for d in atlas.datapoints if d.ligand_restype == ligres.restype]
            size_factor = len(relevant_datapoints)
            for datapoint in relevant_datapoints:
                atoms: Dict[str, Vector] = {}
                for atom, int_coor in datapoint.binder_atoms.items():
                    ext_coor = icoor.internal_to_external(int_coor, False)
                    atoms[atom] = ext_coor
                ares = AlignedResidue(atoms, datapoint.binder_restype, lresid, size_factor)
                aligned_residues.append(ares)
    return LigandBinderComplex(ligand, binder, aligned_residues)
