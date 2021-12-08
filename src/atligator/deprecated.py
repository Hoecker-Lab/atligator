from typing import List

from Bio.PDB import PDBParser, Structure
from Bio.SeqUtils import seq1


def get_aa_from_resnum(structure: str, target: List[int] or int) -> List[str] or str:
    """
    Returns the residue name (one letter code) of a target residue (id or list of ids).
    The name is read out of the structure stored in structure parameter path.
    :param structure: path of structure to use as search space
    :param target: int or list of ints to check for residue name
    :return: residue name or list of residue names corr. to target
    """
    parser = PDBParser(QUIET=True)
    mol: Structure = parser.get_structure(structure, structure)
    if type(target) == int:
        for residue in mol.get_residues():
            if residue.get_id()[1] == target:
                return seq1(residue.get_resname())
        return None
    else:
        resnames = []
        for res_id in target:
            for residue in mol.get_residues():
                if residue.get_id()[1] == res_id:
                    resnames.append(seq1(residue.get_resname()))
    return resnames
