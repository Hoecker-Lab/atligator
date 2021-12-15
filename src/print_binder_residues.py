#!/usr/bin/env python3

"""Utility script that processes a PDB file and identifies binder residues that interact with residues of a ligand,
identified by a chain id.

Example:
./print_binder_residues.py 5AEI.pdb D --radius 3.5 --include_hydrogens

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-03-22
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from operator import itemgetter
from typing import List

from Bio.PDB.Chain import Chain
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from atligator.pdb_util import get_chain_by_id, get_residues_within_radius, is_amino_acid_residue

# Main script. It identifies the chain of the PDB filed specified in the argument list, finds interacting residues of
# non-ligand chains that lie within the specified radius, collects data, and prints them out in a sorted way.
if __name__ == "__main__":
    # noinspection PyTypeChecker
    ap = ArgumentParser(description="processes a specified chain of a PDB file and prints out all residues of "
                                    "different chains that lie within a given interaction radius.",
                        formatter_class=ArgumentDefaultsHelpFormatter,
                        epilog="the output format is: [chain_id] [residue-type][residue_id]")
    ap.add_argument("pdb_file", type=str, help="the PDB file to analyze.")
    ap.add_argument("chain", type=str, help="the relevant chain of the PDB file.")
    ap.add_argument("-r", "--ir_default", type=float, default=4.0,
                    help="pairs of residue atoms within this radius are considered as interacting.")
    ap.add_argument("-b", "--ir_hbond", type=float, default=6.0,
                    help="interaction radius for hydrogen bonded atom pairs.")
    ap.add_argument("-a", "--ir_aromatic", type=float, default=6.0,
                    help="interaction radius for aromatic interaction between atoms.")
    ap.add_argument("-i", "--ir_ionic", type=float, default=8.0,
                    help="interaction radius for ionic binding between positive and negative sidechain atoms.")
    ap.add_argument("-y", "--include_hydrogens", action="store_true",
                    help="whether to include hydrogen atoms. Caveat: many structures do not contain hydrogens.")
    ap.add_argument("-v", "--include_alternative_models", action="store_true",
                    help="whether to include data from models different from the first model of the structure "
                         "(i.e., whether not to ignore alternative NMR conformations)")
    args = ap.parse_args()
    parser = PDBParser(QUIET=True)
    residues = set()
    structure: Structure = parser.get_structure(args.pdb_file, args.pdb_file)
    for model in structure:
        ligand: Chain = get_chain_by_id(model, args.chain)
        for ligres in ligand.get_residues():
            for residue in get_residues_within_radius(model, ligand, ligres, args.ir_default,
                                                      args.ir_hbond, args.ir_aromatic, args.ir_ionic,
                                                      include_hydrogens=args.include_hydrogens):
                if is_amino_acid_residue(residue):
                    residues.add((residue.get_parent().get_id(), residue.get_resname(), int(residue.get_id()[1])))
        if not args.include_alternative_models:
            break
    sorted_residues: List[Residue] = sorted(residues, key=itemgetter(0, 2))
    for entry in sorted_residues:
        print(f"{entry[0]} {entry[1]} {entry[2]}")
