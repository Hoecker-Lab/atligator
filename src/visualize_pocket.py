#!/usr/bin/env python3

"""visualizes the preferred binder residue locations and orientations of a specific ligand residue type based on an
existing atlas file.

Example:
./visualize_pocket.py my.pocket -l ARG -b GLU

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-03-02
"""
import logging
import pickle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import Counter
from typing import List, Dict

from atligator.pdb_util import aa_1to3_conv
from atligator.pocket_miner import Pocket
from atligator.visualization import visualize_pocket_atlas

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# Main script. It processes the atlas file specified in the first argument and visualizes the residues associated with
# the ligand residue types provided in the second argument
if __name__ == "__main__":
    # noinspection PyTypeChecker
    ap = ArgumentParser(description="visualizes the residues in a single pocket.",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("pockets_file", type=str, help="an .pockets file that contains pockets to visualize.")
    ap.add_argument("-l", "--ligand_restype", type=str,
                    help="This defines the ligand residue type. Use three-letter codes, e.g., ARG.")
    ap.add_argument("-p", "--pocket_itemset", type=str,
                    help="Defines the itemset of binder residues of a certain pocket, e.g. TS for a pocket of one Ser "
                         "and one Thr.")
    ap.add_argument("-m", "--method", type=str, default="plotly",
                    help="Atligator offers two plotting methods: 'plotly' or 'matplotlib'")
    args = ap.parse_args()
    with open(args.pockets_file, 'rb') as fo:
        pockets: Dict[str, List[Pocket]] = pickle.load(fo)

    try:
        ligand_pockets = pockets[args.ligand_restype]
    except KeyError as e:
        logger.error(f"The ligand_restype {args.ligand_restype} can not be found in {list(pockets.keys())}")
        raise e

    try:
        desired_itemset = Counter(aa_1to3_conv[res1.upper()] for res1 in args.pocket_itemset)
    except KeyError as e:
        logger.error(f"The pocket_itemset {args.pocket_itemset} includes non-canonical amino acid characters.")
        raise e

    try:
        pocket = [po for po in ligand_pockets if po.itemset.items == desired_itemset][0]
    except IndexError as e:
        logger.error(f"The desired pocket can not be found in these ligand pockets:" + str(ligand_pockets))
        raise e

    visualize_pocket_atlas(pocket,
                           method=args.method)
