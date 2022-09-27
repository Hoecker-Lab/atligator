#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Goes through a list of pdbs and creates new pdbs for each chain longer than min_binder_len residues.
This also contains all other chains within max_distance distance which are parts of segments of >= min_ligand_len
residues.
By default the processed files are written into the current directory. Thus, it is suggested to name an output_path

Should be performed before select_by_sec_struc.py

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-07
"""
import logging
import pathlib
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from glob import iglob
from typing import Dict

from atligator.chain_processing import MultiProcessChainProcessing

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def work(args) -> Dict:
    logger.info(f"Started preprocessing of given pdb structures.")
    start_time = time.perf_counter()
    if args.output_path is None:
        output_path = "./"
    else:
        output_path = args.output_path
    pdbs = [x for x in iglob(args.pdb_files) if ".pdb" in x]
    process = MultiProcessChainProcessing(pdbs=pdbs, output_path=output_path, min_binder_len=args.min_binder_len,
                                          min_ligand_len=args.min_ligand_len, max_distance=args.max_distance,
                                          n_workers=args.n_workers, progress=args.progress_bar)
    objects_dict = process.main()

    with open(pathlib.Path(output_path) / "processed_objects.txt", "w") as objects_file:
        objects_file.write("# The following input files have been processed.\n"
                           f"# Parameters: {args.min_binder_len=} {args.min_ligand_len=} {args.max_distance=}\n")
        for in_file, out_list in objects_dict.items():
            objects_file.write("IN " + in_file + "\n")
            for out_file in out_list:
                objects_file.write(" OUT " + out_file + "\n")

    end_time = time.perf_counter()
    logger.info(f"CPU time wasted for chain processing: {int((end_time - start_time) / 60)} min "
                f"{int((end_time - start_time) % 60)} sec.")
    return objects_dict


if __name__ == "__main__":
    # noinspection PyTypeChecker
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    gr = '\033[37m'
    wt = '\033[0m'
    parser.add_argument("pdb_files", type=str, help=gr + "A pathname of one or more PDB files to "
                                                         "process. Example: './*.pdb'" + wt)
    parser.add_argument("-d", "--max_distance", type=float, default=5.0, help=(gr + "Maximum distance between ligand "
                                                                                    "and binder in angstroms. This will"
                                                                                    " be the maximum distance usable "
                                                                                    "for subsequent atlas generation."
                                                                               + wt))
    parser.add_argument("-o", "--output_path", type=str, default=None, help=(gr + "Output directory. Default is"
                                                                                  "'./'" + wt))
    parser.add_argument("-mb", "--min_binder_len", type=int, default=60, help=(gr + "Threshold which chains to pick as "
                                                                                    "host chains" + wt))
    parser.add_argument("-ml", "--min_ligand_len", type=int, default=3, help=(gr + "Threshold which segments to leave "
                                                                                   "in new pdb" + wt))
    parser.add_argument("-pro", "--progress_bar", action='store_true',
                        help=gr + "Enables a progress bar during processing" + wt)
    parser.add_argument("-w", "--n_workers", type=int, default=-1,
                        help=(gr + "Number of tasks to perform in parallel. -1 corresponds to number of cpus -2." + wt))
    arguments = parser.parse_args()

    work(arguments)
