#!/usr/bin/env python3

"""Performs a naive grafting to a given scaffold structure based on a pocket collection obtained from pocket mining.
Every ligand residue is superimposed with suitable pockets from the collection in order to find the best matching
residue type for the wildcard residues defined in the descriptor. The script also provides a permutation facility, such
that two groups of ligand residues can be assigned wildcards 'XXX' and 'YYY', each of which are replaced by the 20
canonical amino acids, such that 400 combinations can be generated at once.

Example:
./apply_simple_graft.py scaffold.pdb descriptor pockets -o graft.pdb --penalty_threshold 20.0

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-06-13
"""

import pickle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from copy import deepcopy
from typing import Dict, List

from Bio.PDB import PDBParser, Structure, PDBIO

from atligator.acomplex import read_complex_descriptor
from atligator.grafting import simple_graft
from atligator.pdb_util import canonical_amino_acids
from atligator.pocket_miner import Pocket

# Main script. Parses the arguments and calls the simple graft algorithm from the library.
if __name__ == "__main__":
    # noinspection PyTypeChecker
    ap = ArgumentParser(description="analyzes an atlas in terms of frequent pocket itemsets and splits the information"
                                    " into disjoint atlases for every pocket identified.",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("scaffold", type=str,
                    help="a PDB scaffold file containing the binder and ligand coordinates")
    ap.add_argument("descriptor", type=str,
                    help="annotates the scaffold file with information about ligand residue types and wildcard residues"
                         " of the binder. The ligand section may contain the placeholders 'XXX' and 'YYY', which are"
                         " replaced by every of the 20 amino acids in a batch run (see options -x and -y)")
    ap.add_argument("pockets", type=str,
                    help="a dumped pockets file created by the pocket miner. These pockets will be taken into account "
                         "for grafting.")
    ap.add_argument("-d", "--distance_factor", type=float, default=2.0,
                    help="weight factor for distance coefficient of the distance function used for identifying " \
                         "graftable pockets. [1/Angstrom]")
    ap.add_argument("-r", "--orient_factor", type=float, default=1.0,
                    help="weight factor for the primary orientation coefficient (C_alpha->C_beta) of the distance "
                         "function. [1/radians]")
    ap.add_argument("-s", "--secor_factor", type=float, default=1.0,
                    help="weight factor for the secondary orientation coefficient (C_alpha->C_O) of the distance "
                         "function. [1/radians]")
    ap.add_argument("-t", "--penalty_threshold", type=float, default=16.0,
                    help="maximum distance between a wildcard and a grafted pocket residue.")
    ap.add_argument("-x", "--XXX", action='store_true',
                    help='whether to introduce a first batch loop by replacing all ligand residues of type XXX with '
                         'every of the 20 amino acids')
    ap.add_argument('-y', "--YYY", action='store_true',
                    help='whether to introduce a second batch loop based on the placeholder YYY. In total, this '
                         'implies 400 runs of the simple graft algorithm.')
    ap.add_argument("-o", "--output_structure", type=str, default='graft.pdb',
                    help="file name of the PDB structure containing the result of this simple graft. If -x and/or -y"
                         " are set, this name must include XXX and/or YYY, which are replaced by the concrete residue"
                         " types in the output")
    args = ap.parse_args()
    descriptor = read_complex_descriptor(args.descriptor)
    with open(args.pockets, 'rb') as fo:
        pockets: Dict[str, List[Pocket]] = pickle.load(fo)
    for aaX in (canonical_amino_acids if args.XXX else ['XXX']):
        descriptor1 = deepcopy(descriptor)
        if args.XXX:
            descriptor1.ligand_info = [(chain_id, residue_id, aaX if new_residue_type == 'XXX' else new_residue_type)
                                       for (chain_id, residue_id, new_residue_type) in descriptor.ligand_info]
        for aaY in (canonical_amino_acids if args.YYY else ['YYY']):
            descriptor2 = deepcopy(descriptor1)
            out_filename = args.output_structure.replace('XXX', aaX).replace('YYY', aaY)
            print(f'---- generating {out_filename} ----')
            if args.YYY:
                descriptor2.ligand_info = [
                    (chain_id, residue_id, aaY if new_residue_type == 'YYY' else new_residue_type)
                    for (chain_id, residue_id, new_residue_type) in descriptor2.ligand_info]
            parser = PDBParser(QUIET=True)
            scaffold_structure: Structure = parser.get_structure(args.scaffold, args.scaffold)
            simple_graft(scaffold_structure, descriptor2, pockets,
                         penalty_threshold=args.penalty_threshold,
                         distance_factor=args.distance_factor,
                         orient_factor=args.orient_factor,
                         secor_factor=args.secor_factor)
            io = PDBIO()
            io.set_structure(scaffold_structure)
            io.save(out_filename)
