#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Selects processed structures based on secondary structure content.

Should be performed after process_chains.py

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-07
"""
import logging
import pathlib
import shutil
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from glob import iglob
from typing import Union, List

from atligator.construct_selection import SelectionProcess

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def work(args, preprocessed_files: Union[List[str], None] = None):
    logger.info(f"Started selecting pdb structures by given parameters.")
    start_time = time.perf_counter()

    if preprocessed_files is None:
        pdbs = [f for f in iglob(args.pdb_files) if ".pdb" in f]
    else:
        pdbs = preprocessed_files
    process = SelectionProcess(pdbs=pdbs)
    if not pathlib.Path(args.output_path).is_dir():
        pathlib.Path(args.output_path).mkdir()
    filtered = process.filter_pdbs(alphatotal=args.alphatotal, betatotal=args.betatotal, alphabeta=args.alphabeta,
                                   betaalpha=args.betaalpha)
    if args.output_path is None:
        args.output_path = "./"

    with open(pathlib.Path(args.output_path) / "filtered_objects.txt", "w") as objects_file:
        objects_file.write("# The following input files have been filtered.\n"
                           f"# Parameters: {args.alphatotal=} {args.betatotal=} {args.alphabeta=} {args.betaalpha=}\n")
        for file in filtered:
            objects_file.write(file + "\n")

            if args.copy_wanted:
                try:
                    shutil.copy(file, pathlib.Path(args.output_path) / pathlib.Path(file).name)
                except PermissionError as e:
                    logger.warning(e)

    end_time = time.perf_counter()
    logger.info(f"CPU time wasted for selection process: {int((end_time - start_time) / 60)} min "
                f"{int((end_time - start_time) % 60)} sec.")


if __name__ == "__main__":
    # noinspection PyTypeChecker
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    gr = '\033[37m'
    wt = '\033[0m'
    parser.add_argument("pdb_files", type=str, help=gr + "A pathname of one or more PDB files to process." + wt)
    parser.add_argument("-at", "--alphatotal", type=float, default=None, help=gr + "Minimum ratio alpha helix vs total")
    parser.add_argument("-bt", "--betatotal", type=float, default=None, help=gr + "Minimum ratio beta sheet vs total")
    parser.add_argument("-ba", "--betaalpha", type=float, default=None, help=gr + "Minimum ratio beta sheet vs "
                                                                                  "alpha helix")
    parser.add_argument("-ab", "--alphabeta", type=float, default=None, help=gr + "Minimum ratio alpha helix vs "
                                                                                  "beta sheet")
    parser.add_argument("-o", "--output_path", type=str, default=None,
                        help=(gr + "Output directory. Default is './'" + wt))
    parser.add_argument("-cp", "--copy_wanted", action='store_true',
                        help=gr + "Disables copying all filtered structure files in additional '/selected' folder" + wt)

    arguments = parser.parse_args()

    work(arguments)
