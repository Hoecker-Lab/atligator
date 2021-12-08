#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Combines process_chains.py and select_by_sec_struc.py in one application.

Parameters can be used like in both scripts separately. The output path of the first part is used as
input path for select_by_sec_struc.py

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-07-16
"""
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from process_chains import work as process
from select_by_sec_struc import work as select

if __name__ == "__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    gr = '\033[37m'
    wt = '\033[0m'
    parser.add_argument("input_dir", type=str, help=gr + "A pathname of a directory of PDB files to "
                                                         "process" + wt)
    parser.add_argument("-d", "--distance", type=float, default=5.0, help=(gr + "Maximum distance between ligand and "
                                                                                "binder in angstroms" + wt))
    parser.add_argument("-lst", "--pdb_list", type=str, default=None, help=gr + "Optional list of pdb files to "
                                                                                "process" + wt)
    parser.add_argument("-op", "--output", type=str, default=None,
                        help=(gr + "Output directory for processing. Default "
                                   "is 'input_dir/proc/" + wt))
    parser.add_argument("-os", "--output_selection", type=str, default="selected",
                        help=(gr + "Output directory. Default is 'input_dir/../selected/', if processing output is "
                                   "'.../proc/' - else it is 'input_dir/selected" + wt))
    parser.add_argument("-th", "--threshold", type=int, default=60, help=(gr + "Threshold which chains to pick as host"
                                                                               " chains" + wt))
    parser.add_argument("-min", "--minimum", type=int, default=3, help=(gr + "Threshold which segments to leave in new"
                                                                             " pdb" + wt))
    parser.add_argument("-pro", "--progress_bar", action='store_true',
                        help=gr + "Enables a progress bar during processing" + wt)
    parser.add_argument("-w", "--n_workers", type=int, default=-1,
                        help=(gr + "Number of tasks to perform in parallel. -1 corresponds to number of cpus -2." + wt))

    # selection:
    parser.add_argument("-at", "--alphatotal", type=float, default=0.55, help=gr + "Minimum ratio alpha helix vs total")
    parser.add_argument("-ba", "--betaalpha", type=float, default=0.1, help=gr + "Maximum ratio beta sheet vs "
                                                                                 "alpha helix")
    parser.add_argument("-cp", "--copy_wanted", action='store_false',
                        help=gr + "Disables copying all wanted object files in additional '/selected' folder" + wt)
    parser.add_argument("-pv", "--pymol_visuals", action='store_true',
                        help=gr + "Enables visualization in pymol" + wt)
    args = parser.parse_args()

    args.input_dir = process(args)
    #args.output = args.output_selection
    #select(args)
