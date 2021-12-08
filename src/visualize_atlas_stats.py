#!/usr/bin/env python3

"""Visualizes statistical information about an atlas as stacked bar plot. This shows the number of datapoints of every
combination of ligand and residue type.

Example:
./visualize_atlas_stats.py my_atlas

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-07-13
"""

import pickle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType

from atligator.visualization import visualize_atlas_stats

# Main script. It opens the atlas file specified in the first argument and calls the atligator library function for
# visualizing its statistics.
if __name__ == "__main__":
    ap = ArgumentParser(description="visualizes statistical information about an atlas.",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("atlas_file", type=FileType('rb'), help="an .atlas file to visualize.")
    ap.add_argument("-i", "--ligand_per_binder", action="store_true",
                    help="If enabled: inverts the representation, i.e., group ligand residues by binder residues")
    args = ap.parse_args()
    with open(args.atlas_file.name, 'rb') as fo:
        atlas = pickle.load(fo)
    visualize_atlas_stats(atlas, args.ligand_per_binder)
