#!/usr/bin/env python3

"""Generates a dataset of positional properties for specific ligand-binder interaction.

Example:
./generate_atlas.py my.atlas structures/*.pdb --minlen 4 --maxlen 30 --ir_default 3.5

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-02-26
"""

import glob
import pickle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from typing import List

from atligator.atlas import generate_atlas

# Main script. Processes the files specified in the argument list (considering also wildcards), identifies atlas
# datapoints in there, combines them into an atlas, and prints the atlas out.
if __name__ == "__main__":
    ap = ArgumentParser(description="creates an atlas of specific ligand-binder interaction properties from a library "
                                    "of PDB files.",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("output_file", type=FileType('wb'),
                    help="path to the binary output file. Typical extension is '.atlas'.")
    ap.add_argument("pdb_files", type=str, nargs="+", help="a list of pathnames or pathname patterns of PDB files")
    ap.add_argument("-n", "--minlen", type=int, default=4,
                    help="the minimum length of a chain to be treated as ligand")
    ap.add_argument("-x", "--maxlen", type=int, default=20,
                    help="the maximum length of a chain to be treated as ligand")
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
    ap.add_argument("-l", "--alternative_lig_aa", type=str, default="",
                    help="If given this alternative name of non-canonical aas will be also searched for as a ligand residue")
    ap.add_argument("-o", "--restrict_to_alternative", action="store_true",
                    help="If True: If alternative Ligand amino acid is given restrict to this type.")
    ap.add_argument("-s", "--allow_self_interactions", action="store_true",
                    help="If True interactions within one chain are taken into account (Slows down calculation heavily!!)")
    ap.add_argument("-sb", "--skip_bb_atoms", action="store_true",
                    help="If True interactions of backbone atoms are not considered (recommended!)")
    ap.add_argument("-v", "--include_alternative_models", action="store_true",
                    help="whether to include data from models different from the first model of the structure "
                         "(i.e., whether not to ignore alternative NMR conformations)")
    ap.add_argument("-k", "--n_workers", type=int, default=-1,
                    help="Number of workers to use for parallel processing. -1 for number of CPU cores")
    args = ap.parse_args()
    filenames: List[str] = []
    for arg in args.pdb_files:
        filenames.extend(glob.glob(arg))

    total = len(filenames)


    def funny_observer(_name: str):

        print(f"Done with {filenames.index(_name) + 1} of {total}\t\t{_name}")


    atlas = generate_atlas(filenames, args.ir_default, args.ir_hbond, args.ir_aromatic, args.ir_ionic,
                           args.minlen, args.maxlen, args.n_workers, args.include_hydrogens,
                           args.include_alternative_models,
                           alternative_lig_aa=args.alternative_lig_aa if args.alternative_lig_aa else None,
                           restrict_to_alternative=args.restrict_to_alternative,
                           allow_self_interactions=args.allow_self_interactions, observer=funny_observer,
                           skip_bb_atoms=args.skip_bb_atoms)
    print("Generated atlas with", len(atlas.datapoints), "Datapoints.")
    with open(args.output_file.name, 'wb') as fo:
        pickle.dump(atlas, fo)
    print("Saved Atlas in", args.output_file.name)
