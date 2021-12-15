# -*- coding: utf-8 -*-
"""
For preprocessing chains (multithreaded).
Each host chain (> threshold!) is processed individually to remove DNA/RNA, H2O and short segments.
All chains other than host chain get segmented and called with individual chain ids. Segments shorter than
minimum get discarded.
One pdb file with this construct is saved.

Update 8-8-18:
Removed htmd dependencies completely and speeded up calculation. a.118.1 (511 structures, 560 MB) takes 3 min 30 sec.
On 4 cores of i7 7700K.


:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-04-05
"""
import logging
import os
import pathlib
# for multithreading without progress bar:
# ProcessPoolExecutor better than Thread..." 50% time reduction (without progress bar option)
from concurrent.futures import ProcessPoolExecutor as ProcessPool
from copy import copy
from multiprocessing import Pool
from typing import Any, Dict, List, Set, Tuple

from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB import Residue
# For progress bar during multiprocessing:
from tqdm import tqdm

from atligator.pdb_util import ChainSelect, NoBinderAtoms, NoLigandAtoms, NoResidueInDistance
from atligator.pdb_util import get_nonbinder_residues_within_radius

logger = logging.getLogger(__name__)

def multiprocess(array, func, n=1) -> Dict:
    """
    Perform function on array with n parallel processes.
    :param array: Array with parameters for function.
    :param func: Function which is executed with array parameters.
    :param n: Number of processes
    :return: list of output of function
    """
    with ProcessPool(max_workers=n) as pool:
        return pool.map(func, array)


def process_threaded_w_prog(array, func, n=1) -> Dict:
    """
    Perform function on array with n parallel processes. Includes a progress bar.
    :param array: Array with parameters for function (file paths!).
    :param func: Function which is executed with array parameters.
    :param n: Number of processes
    :return: dictionary of output of function
    """

    sizecounter: int = 0

    element_size = {}
    for file in array:
        size = os.stat(file).st_size
        sizecounter += size
        element_size[file] = size

    if sizecounter >= 1024:
        print_size = "{0:.2f}".format(sizecounter / 1024) + " kB"
        if sizecounter >= 1024 * 1024:
            print_size = "{0:.2f}".format(sizecounter / 1024 / 1024) + " MB"
    else:
        print_size = f"{sizecounter} B"

    logger.info(f"Found {len(array)} files with a total size of {print_size}")
    return imap_unordered(sizecounter, array, element_size, func, n)


def wrap_my_func(args: Tuple[int, callable]) -> Tuple[int, Any]:
    """
    Wraps a function (process_pdb) and returns result
    :param args: Tuple of an int to be processed by function and function.
    :return: Tuple of int parameter and result of function called with int.
    """
    res = args[1](args[0])
    return args[0], res


def imap_unordered(sizecounter, array, elementsize, func, n):
    with tqdm(total=sizecounter, unit='B', unit_scale=True) as pbar:
        output_dict = {}
        pool = Pool(processes=n)
        arg = ((i, func) for i in array)
        # map / imap is only executed when return is used. same for apply_async: use .get()!
        for i, res in pool.imap_unordered(wrap_my_func, arg):
            output_dict[i] = res
            pbar.update(elementsize[i])
        pool.close()
        pool.join()
        return output_dict


def process_pdb(pdb, max_distance, min_ligand_len, min_binder_len, output_path, quiet: bool = True) -> List[str]:
    """
    Processes a pdb file.
    First checks chain length of all chains. If chain is more than min_binder_len, it will be processed
    furtheron as a host chain.
    For each host chain all other chains are removed, if not in distance max_distance of host.
    All non-amino acid residue also get removed.
    The remaining chain other than the host chain get segmented. Segments are removed if they don't consist
    of at least 3 residues (min_ligand_len). If theres no other chain than the host, it gets ignored.
    Molecules with host chain and segments are saved as new pdb files with
    a letter for the chain name.
    This name is also returned.
    :param pdb: The pdb file (OpenFileHandle) or pdb file name (relative or absolute path) to open
    :param max_distance: Maximum distance of host and ligand chain
    :param min_ligand_len: Number of residues a ligand chain must have at least to be processed
    :param min_binder_len: Number of residues a chain must have at least to act as a host chain
    :param output_path: The path where to save the processed files
    :param quiet: If true, all prints are suppressed
    :return: List of the names of the newly created files.
    """
    if not isinstance(pdb, str):
        pdb_name = "structure"
    else:
        pdb_name = pdb.replace(".pdb", "")

    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    mol = parser.get_structure(id=pdb_name, file=pdb)

    return_list = []

    chain_counter = 0
    for models in mol:
        for chain in models:
            if len(chain) >= min_binder_len:
                # Find all the residues in distance of the binder chain
                in_distance: Set[Residue] = set()
                try:
                    for res_in_dist in (get_nonbinder_residues_within_radius(models, chain,
                                                                             ir_max=max_distance,
                                                                             ligand_min_len=min_ligand_len,
                                                                             include_hydrogens=False)):
                        in_distance.add(res_in_dist)
                except (NoBinderAtoms, NoLigandAtoms, NoResidueInDistance):
                    pass  # in_distance is still empty

                if in_distance:
                    molch = copy(mol)
                    io = PDBIO()
                    io.set_structure(molch)
                    io.save(str(pathlib.Path(output_path) / pathlib.Path(pdb_name).name) + f"{chain.get_id()}.pdb",
                            ChainSelect(in_distance, chain.get_id()))
                    return_list.append(
                        str(pathlib.Path(output_path) / pathlib.Path(pdb_name).name) + f"{chain.get_id()}.pdb"
                    )
                    chain_counter += 1
                if not quiet:
                    print(f"step {chain_counter} out of {len(list(mol.get_chains()))} for {pdb_name}, ")
    return return_list


class MultiProcessChainProcessing:
    """
    Performs the ALARMS processing step for a set of pdb structures.
    For multiprocessing it makes use of either
    - multiprocess() for processing without a progress bar
    or
    - process_threaded_w_prog() for processing with a progress bar
    """

    def __init__(self, pdbs: List[str], output_path: str, min_binder_len: int, min_ligand_len: int, max_distance: float,
                 n_workers: int = -1, progress: bool = True, quiet: bool = False):
        """
        :param pdbs: A list of pdb structures to process
        :param output_path: Paths to crucial directories (input, output)
        :param min_binder_len: Number of residues a chain must have at least to act as a host chain
        :param min_ligand_len: Number of residues a ligand chain must have at least to be processed
        :param max_distance: Maximum distance of host and ligand chain
        :param n_workers: Number of processes
        :param progress: If true, a progress bar is displayed (tqdm)
        :param quiet: If true, all prints are suppressed
        """
        self.pdbs = pdbs
        self.output_path = output_path
        if not pathlib.Path(output_path).is_dir():
            pathlib.Path(output_path).mkdir()
        self.min_binder_len = min_binder_len
        self.min_ligand_len = min_ligand_len
        self.max_distance = max_distance
        if n_workers < 0:
            self.n_workers = os.cpu_count() - 2
        else:
            self.n_workers = n_workers
        self.progress = progress
        self.QUIET = quiet if not progress else True

    def process_pdb(self, pdb):
        return process_pdb(pdb=pdb, max_distance=self.max_distance, min_ligand_len=self.min_ligand_len,
                           min_binder_len=self.min_binder_len, output_path=self.output_path, quiet=self.QUIET)

    def main(self) -> Dict:
        """
        creates set of processed pdbs to save in objects.txt file.
        :return : A set of pdb names with their corresponding main chain, e.g. 2GDNA with '2GDN' as pdb identifier and
        A as the main chain in this object.
        """
        if self.progress:
            objects_dict = process_threaded_w_prog(self.pdbs, self.process_pdb, self.n_workers)
        else:
            objects_dict = {k: v for k, v in zip(self.pdbs, multiprocess(self.pdbs, self.process_pdb, self.n_workers))}
        return objects_dict
