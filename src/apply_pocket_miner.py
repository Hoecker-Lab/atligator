#!/usr/bin/env python3

"""Identifies frequent binding pockets in an atlas. The result, a collection of several atlases per ligand residue
type, is exported into a specified target directory. The data points contained in the individual atlases are
clustered based on an extended k-means clustering algorithm. Optionally, the script also generates structures and/or
HTML pages summarizing the calculated information (and structures).

Example:
./apply_pocket_miner.py my.atlas --confidence_threshold 0.1 --output_directory pockets -c -p -cp -y -a 8

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-05-18
"""

import pickle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from pickle import load, dump

from Bio.PDB.PDBIO import PDBIO
from Bio.SeqUtils import seq1

from atligator.html_util import generate_html_pocket_index, generate_html_pocket_description
from atligator.pocket_miner import mine_pockets

# Main script. Processes the atlas, analyzes it in terms of pocket itemsets, identifies frequent pocket itemsets,
# and uses this information to filter the atlas by different pockets, such that one individual atlas per pocket
# detected is obtained
if __name__ == "__main__":
    ap = ArgumentParser(description="analyzes an atlas in terms of frequent pocket itemsets and splits the information"
                                    " into disjoint atlases for every pocket identified.",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("atlas", type=FileType('r'), help="an atlas to serve as input for the analysis")
    ap.add_argument("-x", "--max_atlases", type=int, default=5,
                    help="maximum number of atlases generated for every ligand residue type")
    ap.add_argument("-m", "--min_cardinality", type=int, default=3,
                    help="minimum cardinality of itemsets, i.e. the number of residues part of a pocket")
    ap.add_argument("-n", "--confidence_threshold", type=float, default=0.05,
                    help="confidence threshold for association rule mining")
    ap.add_argument("-t", "--support_threshold", type=float, default=0.01,
                    help="support threshold for association rule mining")
    ap.add_argument("-d", "--distance_factor", type=float, default=1.0,
                    help="weight factor for distance coefficient of the distance function used for clustering."
                         " [1/Angstrom]")
    ap.add_argument("-r", "--orient_factor", type=float, default=2.0,
                    help="weight factor for the primary orientation coefficient (C_alpha->C_beta) of the distance "
                         "function used for clustering. [1/radians]")
    ap.add_argument("-s", "--secor_factor", type=float, default=1.0,
                    help="weight factor for the secondary orientation coefficient (C_alpha->C_O) of the distance "
                         "function used for clustering. [1/radians]")
    ap.add_argument("-v", "--variance_threshold", type=float, default=5.0,
                    help="the residues contained by every cluster are allowed to have at most this variance; if this "
                         "is exceeded, the cluster is reduced in an iterative post-processing step")
    ap.add_argument("-o", "--output_folder", type=str, default='pockets',
                    help="the created pocket atlases will be stored in this folder. This directory must be existant "
                         "and the path must not include a trailing '/'.")
    ap.add_argument("-f", "--pockets_filename", type=str, default='pockets.pockets',
                    help="name of the to be dumped pocket file inside the output folder")
    ap.add_argument("-a", "--save_atlases", action="store_true",
                    help="whether to save the individual clusters for ligand residues as separate atlas files "
                         "(suffix .atlas)")
    ap.add_argument("-ac", "--save_atlas_centroids", action="store_true",
                    help="whether to save the centroids of the clusters into a separate atlas file (suffix "
                         "_clustered.atlas).")
    ap.add_argument("-p", "--save_pdb", action="store_true",
                    help="whether to export the clusters into PDB structures, where the residues part of the clusters "
                         "are represented in a superimposed way.")
    ap.add_argument("-pc", "--save_pdb_centroids", action="store_true",
                    help="whether to save the centroids of the clusters into a separate PDB structure file (suffix "
                         "_clustered.pdb).")
    ap.add_argument("-y", "--generate_html", action="store_true",
                    help="whether to generate a HTML based static website that shows information about the generated "
                         "atlas and visualizes saved PDB structures using the NGL.js renderer.")
    ap.add_argument("-k", "--n_workers", type=int, default=-1,
                    help="Number of workers to use for parallel processing. -1 for number of CPU cores")
    args = ap.parse_args()
    with open(args.atlas.name, 'rb') as fo:
        atlas = load(fo)
    pockets = mine_pockets(atlas,
                           max_atlases_per_ligand_restype=args.max_atlases,
                           min_itemset_cardinality=args.min_cardinality,
                           confidence_threshold=args.confidence_threshold,
                           support_threshold=args.support_threshold,
                           distance_factor=args.distance_factor,
                           orient_factor=args.orient_factor,
                           secor_factor=args.secor_factor,
                           variance_threshold=args.variance_threshold,
                           n_workers=args.n_workers)
    with open(args.output_folder + '/' + args.pockets_filename, 'wb') as fo:
        pickle.dump(pockets, fo)
    for lig_restype, pocketlist in pockets.items():
        for pocket in pocketlist:
            filename = f"{args.output_folder}/{seq1(lig_restype)}_{pocket.pocket_id()}_" \
                       f"{int(pocket.support * 100.0)}.atlas"
            pocket_atlas = pocket.to_atlas()
            if args.save_atlases:
                with open(filename, 'wb') as fo:
                    dump(pocket_atlas, fo)
            if args.save_atlas_centroids:
                clustered_filename = filename.replace(".atlas", "_clustered.atlas")
                clustered_atlas = pocket.to_clustered_atlas(args.distance_factor, args.orient_factor, args.secor_factor)
                with open(clustered_filename, 'wb') as fo:
                    dump(clustered_atlas, fo)
            if args.save_pdb:
                structure_filename = filename.replace(".atlas", ".pdb")
                pocket_structure = pocket.to_pdb_structure()
                io = PDBIO()
                io.set_structure(pocket_structure)
                io.save(structure_filename)
            if args.save_pdb_centroids:
                clustered_structure_filename = filename.replace(".atlas", "_clustered.pdb")
                clustered_pocket_structure = pocket.to_clustered_pdb_structure(args.distance_factor, args.orient_factor,
                                                                               args.secor_factor)
                io = PDBIO()
                io.set_structure(clustered_pocket_structure)
                io.save(clustered_structure_filename)
        if args.generate_html:
            generate_html_pocket_description(args.output_folder, lig_restype, pocketlist, args.save_pdb,
                                             args.save_pdb_centroids)
    if args.generate_html:
        generate_html_pocket_index(args.output_folder, atlas, pockets)
