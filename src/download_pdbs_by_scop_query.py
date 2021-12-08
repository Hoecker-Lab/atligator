#!/usr/bin/env python3

"""Finds all PDB files that match a given SCOP(e) query using a local copy of the SCOP(e) database, downloads the
structures from a local copy of the PDB database, and stores the files in PDB format under the given destination
directory.

Example:
./download_pdbs_by_scop_query a.118.1 -s /agh/db/scop/2.07/dir.cla.scope.2.07-stable.txt -p /agh/db/pdb

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-06-11
"""
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType

# Main script. Parses the arguments and passes them to the library function.
from atligator.pdb_util import download_pdbs_by_scop_query

if __name__ == "__main__":
    ap = ArgumentParser(description="Extracts PDB structures matching a SCOP query to a specified destination folder.",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("query", type=str, help="A SCOPe fold, family, or superfamily. E.g., 'a.118.1'")
    ap.add_argument("-s", "--scop_dir_cla_file", type=FileType('r'),
                    default='/agh/db/scop/2.07/dir.cla.scope.2.07-stable.txt',
                    help="part of the SCOP(e) database that contains the classification")
    ap.add_argument("-p", "--pdb_database", type=str, default='/agh/db/pdb',
                    help="directory that contains a copy of the PDB database. Must contain a folder 'data'.")
    ap.add_argument("-o", "--output_folder", type=str, default='structures',
                    help="the extracted structures will be stored in this folder. This directory must be existant "
                         "and the path must not include a trailing '/'.")
    args = ap.parse_args()
    download_pdbs_by_scop_query(query=args.query,
                                dest_dir=args.output_folder,
                                pdb_db_path=args.pdb_database,
                                scop_dir_cla_file=args.scop_dir_cla_file.name)
