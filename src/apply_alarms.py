#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Combines process_chains.py and select_by_sec_struc.py in one application.

Parameters can be used like in both scripts separately. The output path of the first part is used as
input path for select_by_sec_struc.py

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-07-16
"""
import logging
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from process_chains import work as process
from select_by_sec_struc import work as select

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":
    # noinspection PyTypeChecker
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    gr = '\033[37m'
    wt = '\033[0m'
    parser.add_argument("pdb_files", type=str, help=gr + "A pathname of one or more PDB files to "
                                                         "process. Example: './*.pdb'" + wt)
    parser.add_argument("-d", "--max_distance", type=float, default=5.0,
                        help=(gr + "Maximum distance between ligand and binder in angstroms. This will be the maximum "
                                   "distance usable for subsequent atlas generation." + wt))
    parser.add_argument("-op", "--output_path_processing", type=str, default=None,
                        help=(gr + "Output directory for processing. Default is './'" + wt))
    parser.add_argument("-mb", "--min_binder_len", type=int, default=60, help=(gr + "Threshold which chains to pick as "
                                                                                    "host chains" + wt))
    parser.add_argument("-ml", "--min_ligand_len", type=int, default=3, help=(gr + "Threshold which segments to leave "
                                                                                   "in new pdb" + wt))
    parser.add_argument("-pro", "--progress_bar", action='store_true',
                        help=gr + "Enables a progress bar during processing" + wt)
    parser.add_argument("-w", "--n_workers", type=int, default=-1,
                        help=(gr + "Number of tasks to perform in parallel. -1 corresponds to number of cpus -2." + wt))
    parser.add_argument("-at", "--alphatotal", type=float, default=None, help=gr + "Minimum ratio alpha helix vs total")
    parser.add_argument("-bt", "--betatotal", type=float, default=None, help=gr + "Minimum ratio beta sheet vs total")
    parser.add_argument("-ba", "--betaalpha", type=float, default=None, help=gr + "Minimum ratio beta sheet vs "
                                                                                  "alpha helix")
    parser.add_argument("-ab", "--alphabeta", type=float, default=None, help=gr + "Minimum ratio alpha helix vs "
                                                                                  "beta sheet")
    parser.add_argument("-os", "--output_path_selection", type=str, default=None,
                        help=(gr + "Output directory for selection. Default is './'" + wt))
    parser.add_argument("-cp", "--copy_wanted", action='store_true',
                        help=gr + "Disables copying all filtered structure files in additional '/selected' folder" + wt)
    args = parser.parse_args()

    args.output_path = args.output_path_processing
    processed_dict = process(args)
    preprocessed_files = []
    for output_files in processed_dict.values():
        preprocessed_files.extend(output_files)
    args.output_path = args.output_path_selection
    select(args, preprocessed_files=preprocessed_files)
