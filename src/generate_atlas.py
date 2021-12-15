#!/usr/bin/env python3

"""Generates a dataset of positional properties for specific ligand-binder interaction.

Example:
./generate_atlas.py my.atlas structures/*.pdb --min_ligand_len 4 --max_ligand_len 30 --ir_default 3.5

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-02-26
"""

import glob
import logging
import pickle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from typing import List

from atligator.atlas import AtlasGeneration

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# Main script. Processes the files specified in the argument list (considering also wildcards), identifies atlas
# datapoints in there, combines them into an atlas, and prints the atlas out.
if __name__ == "__main__":
    # noinspection PyTypeChecker
    ap = ArgumentParser(description="creates an atlas of specific ligand-binder interaction matches from a collection "
                                    "of PDB files.",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("output_file", type=str,
                    help="Path to the binary output file. Typical extension is '.atlas'.")
    ap.add_argument("pdb_files", type=str, help="A pathname or pathname pattern of PDB files")
    ap.add_argument("-n", "--min_ligand_len", type=int, default=3,
                    help="The minimum length of a chain to be treated as ligand")
    ap.add_argument("-x", "--max_ligand_len", type=int, default=20,
                    help="The maximum length of a chain to be treated as ligand")
    ap.add_argument("-dr", "--ir_default", type=float, default=4.0,
                    help="Pairs of residue atoms within this radius are considered as interacting.")
    ap.add_argument("-hr", "--ir_hbond", type=float, default=6.0,
                    help="Interaction radius for potentially hydrogen bonded atom pairs.")
    ap.add_argument("-ar", "--ir_aromatic", type=float, default=6.0,
                    help="Interaction radius for potentially aromatic interactions between atoms.")
    ap.add_argument("-ir", "--ir_ionic", type=float, default=8.0,
                    help="Interaction radius for potential charged interactions between positively and negatively "
                         "charged sidechain atoms.")
    ap.add_argument("-y", "--include_hydrogens", action="store_true",
                    help="Whether to include hydrogen atoms. Caveat: Many input structures do not contain hydrogens.")
    ap.add_argument("-a", "--alternative_lig_aa", type=str, default="",
                    help="If given this alternative name of non-canonical amino acids will be also searched for as a "
                         "ligand residue")
    ap.add_argument("-r", "--restrict_to_alternative", action="store_true",
                    help="If True: If alternative Ligand amino acid is given restrict to this type.")
    ap.add_argument("-s", "--allow_self_interactions", action="store_true",
                    help="If True interactions within one chain are taken into account (slows down calculation "
                         "heavily!!)")
    ap.add_argument("-b", "--skip_bb_atoms", action="store_true",
                    help="If True interactions of backbone atoms are not considered (recommended!)")
    ap.add_argument("-m", "--include_alternative_models", action="store_true",
                    help="Whether to include data from models different from the first model of the structure "
                         "(i.e., whether not to ignore alternative NMR structure conformations)")
    ap.add_argument("-w", "--n_workers", type=int, default=-1,
                    help="Number of workers to use for parallel processing. -1 for number of CPU cores")
    args = ap.parse_args()

    logger.info(f"Starting generation of atlas.")
    filenames: List[str] = glob.glob(args.pdb_files)

    total = len(filenames)

    generator = AtlasGeneration(filenames)
    atlas = generator.generate_atlas(ir_default=args.ir_default,
                                     ir_hbond=args.ir_hbond,
                                     ir_aromatic=args.ir_aromatic,
                                     ir_ionic=args.ir_ionic,
                                     min_ligand_len=args.min_ligand_len,
                                     max_ligand_len=args.max_ligand_len,
                                     n_workers=args.n_workers,
                                     include_hydrogens=args.include_hydrogens,
                                     include_alternative_models=args.include_alternative_models,
                                     alternative_lig_aa=args.alternative_lig_aa if args.alternative_lig_aa else None,
                                     restrict_to_alternative=args.restrict_to_alternative,
                                     allow_self_interactions=args.allow_self_interactions,
                                     skip_bb_atoms=args.skip_bb_atoms)
    logger.info(f"Generated atlas with {len(atlas.datapoints)} Datapoints.")
    with open(args.output_file, 'wb') as fo:
        pickle.dump(atlas, fo)
    logger.info(f"Saved Atlas in {args.output_file}")
