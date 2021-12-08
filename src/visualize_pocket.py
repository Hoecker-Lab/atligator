#!/usr/bin/env python3

"""visualizes the preferred binder residue locations and orientations of a specific ligand residue type based on an
existing atlas file.

Example:
./visualize_pocket.py my.pocket -l ARG -b GLU

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-03-02
"""

import pickle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType

from atligator.visualization import visualize_atlas

# Main script. It processes the atlas file specified in the first argument and visualizes the residues associated with
# the ligand residue types provided in the second argument
if __name__ == "__main__":
    ap = ArgumentParser(description="visualizes the residues aligned with a specific ligand type in an atlas.",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("pockets_file", type=FileType('rb'), help="an .pockets file that contains pockets to visualize.")
    ap.add_argument("-l", "--ligand_restype", type=str,
                    help="if set, this restricts the ligand residue type whose aligned binder residues to show "
                         "to a specfic type. Use three-letter codes, e.g., ARG.")
    ap.add_argument("-p", "--pocket_itemset", type=str,
                    help="if set, this restricts the ligand residue type to show to a specfic type. Use three-letter "
                         "codes, e.g., ARG.")
    ap.add_argument("-m", "--method", type=str, default="plotly",
                    help="Atligator offers two plotting methods: 'plotly' or 'matplotlib'")
    ap.add_argument("-2", "--draw_secor", action="store_true",
                    help="Whether to draw secondary orientation vectors of binder residues. Only used with "
                         "matplotlib method.")
    args = ap.parse_args()
    with open(args.atlas_file.name, 'rb') as fo:
        atlas = pickle.load(fo)
    visualize_atlas(atlas,
                    method=args.method,
                    ligand_restype=args.ligand_restype,
                    binder_restype=args.binder_restype,
                    draw_secor=args.draw_secor)
