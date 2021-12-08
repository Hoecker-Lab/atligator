#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Selects processed structures based on secondary structure content.

Should be performed after process_chains.py

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-07
"""
import pathlib
import shutil
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from glob import iglob

from alarms.construct_selection import SelectionProcess


def work(args):
    start_time = time.perf_counter()

    pdbs = [f for f in iglob(args.pdb_files) if ".pdb" in f]
    process = SelectionProcess(pdbs=pdbs)
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
                    print(e)

    end_time = time.perf_counter()
    print("CPU time wasted for selection process:", int((end_time - start_time) / 60), "min",
          int((end_time - start_time) % 60), "sec.")


if __name__ == "__main__":
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
