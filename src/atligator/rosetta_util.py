"""This module abstracts from the Rosetta program suite by wrapping the interface of the relax application.

:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-04-04
"""

import os
import shutil
from glob import glob
from subprocess import run, CalledProcessError
from typing import List, Any

from Bio.PDB.PDBParser import PDBParser
from numpy import std

from atligator.acomplex import ComplexDescriptor
from atligator.base_scorer import BaseScorer, ScoreResult
from atligator.pdb_util import SecondaryScores, inspect_pdb, get_ligand_positions, get_binder_positions


def run_command(rosetta_path: str, exe_file: str, db_path: str, input_file: str, output_dir: str, score_file: str,
                extra_rotamers: bool, mode: str, nstruct: int, script: str or None, stdout: Any,
                resfile: str = "./resfile.res") -> None:
    command = f"{rosetta_path}{os.sep}{exe_file} -database {db_path}" + \
              f" -in:file:s {input_file} -out:path:all {output_dir}" + \
              f" -out:path:score {output_dir} -out:file:scorefile {score_file}" + \
              (" -ex1 -ex2" if extra_rotamers else "")
    if mode == "fixbb":
        command += f" -resfile {resfile} -minimize_sidechains"
    elif mode == "score":
        command += f" -out:nstruct {nstruct}" + (" -relax:fast" if script is None else f" -relax:script {script}")
    try:
        run(command, shell=True, check=True, stdout=stdout)
    except CalledProcessError as err:
        print(f"Unexpected {err=}, {type(err)=}")
        raise


class RosettaScorer(BaseScorer):
    """This class encapsulates Rosetta's relax functionality and provides a minimal scripting interface."""

    def __init__(self, rosetta_path: str = '$ROSETTA', relax_executable: str = 'relax.default.linuxgccrelease',
                 rosetta_db_path: str = '$ROSETTA_DB', mode: str = 'relax', relax_script: str = None,
                 extra_rotamers: bool = False, nstruct: int = 1, ligand_factor: float = 20, lin_factor: float = 30,
                 pen_factor: float = 10.0, rmsd_slope: float = 5, rmsd_offset: float = 5, mutant_dir: str = './',
                 binding_site_factor: float = None):
        """Creates a new instance based on information about the Rosetta executable
        :param rosetta_path: the path to the folder where Rosetta's executable binaries are located
        :param relax_executable: the complete name of the relax executable file (including extensions)
        :param rosetta_db_path: the path to the Rosetta databases
        :param mode: Rosetta structure relaxation mode (one of 'score', 'relax')
        :param relax_script: path to a Rosetta relax script, or None if default protocol shall be kept
        :param extra_rotamers: if true, extra rotamers (ex1, ex2) are used
        :param nstruct: the number of structures to produce in every relax run
        :param ligand_factor: the share of ligand energy in the calculation of total energy
        :param lin_factor: linearity impact in total energy calculation
        :param pen_factor: impact of RMSD penalty
        :param rmsd_slope: The slope of the influence of ligand rmsd after the offset
        :param rmsd_offset: An offset value where rmsd penalty will start
        :param mutant_dir: directory where mutant files are stored - relative to current directory.
        :param binding_site_factor: the share of binding site energy in the calculation of total energy
        """
        self.rosetta_path = rosetta_path
        self.relax_executable = relax_executable
        self.rosetta_db_path = rosetta_db_path
        self.mode = mode
        self.relax_script = relax_script
        self.extra_rotamers = extra_rotamers
        self.nstruct = nstruct
        self.ligand_factor = ligand_factor
        self.binding_site_factor = binding_site_factor if binding_site_factor is not None else ligand_factor
        self.lin_factor = lin_factor
        self.pen_factor = pen_factor
        self.rmsd_slope = rmsd_slope
        self.rmsd_offset = rmsd_offset
        self.mutant_dir = mutant_dir if mutant_dir.endswith('/') else mutant_dir + '/'

    def read_score(self, scorefile_path: str, descriptor: ComplexDescriptor, div_stddev: bool = False) -> ScoreResult:
        """Creates a ScoreResult instance from a score file created by Rosetta.
        :param scorefile_path: the name of the scorefile
        :param descriptor: forwarded complex descriptor
        :param div_stddev: If enabled ligand score will be divided by the standard deviation of indiv. residue scores.
        This is not used atm. (see score() return object).
        :return: the ScoreResult runtime object containing both parsed Rosetta scores and secondary scores
        """
        with open(scorefile_path, 'r') as sf:
            line = ""
            while not line.startswith("SCORE:"):
                line = sf.readline()  # ignore headers (if present)
            best_score = float("inf")
            best_binder_score = float("inf")
            best_ligand_score = float("inf")
            best_binding_site_score = float("inf")
            best_lin_score = float("inf")
            best_rmsd = float("inf")
            best_id = ""
            line = sf.readline()
            while not line == "":
                values = line.split()
                total_score = float(values[1])
                mutant_id = self.mutant_dir + values[-1]
                lig_score: float = 0.0
                bin_site_score: float = 0.0
                std_dev_lst: List[float] = []
                parser = PDBParser(QUIET=True)
                mutant_structure = mutant_id + '.pdb'
                biopdb_structure = parser.get_structure(mutant_structure, mutant_structure)
                ligand_positions = get_ligand_positions(biopdb_structure,
                                                        descriptor.ligand_info)
                binding_site_positions = get_binder_positions(biopdb_structure,
                                                              descriptor.binder_info)
                with open(mutant_structure, 'r') as rf:
                    for line in rf.readlines():
                        if any("_" + str(pos) in line for pos in ligand_positions):
                            std_dev_lst.append(float(line.split()[-1]))
                            lig_score += float(line.split()[-1])
                        elif any("_" + str(pos) in line for pos in binding_site_positions):
                            bin_site_score += float(line.split()[-1])
                secondary_scores: SecondaryScores = inspect_pdb(mutant_id, descriptor.ligand_info)
                lin_score = secondary_scores.linearity_p_res
                if div_stddev:
                    lig_score = lig_score / std(std_dev_lst)
                rmsd = secondary_scores.rmsd
                pen = (rmsd - self.rmsd_offset) ** self.rmsd_slope if rmsd >= self.rmsd_offset else 0.0
                # binding site score must be normalized to fit influence of ligand score
                norm_bin_site_score = bin_site_score * len(descriptor.ligand_info) / len(descriptor.binder_info)
                fitness_score: float = total_score + self.ligand_factor * lig_score + self.lin_factor * lin_score + \
                                       self.pen_factor * pen + self.ligand_factor * norm_bin_site_score
                if fitness_score < best_score:
                    best_score = fitness_score
                    best_binder_score = total_score - lig_score
                    best_ligand_score = lig_score
                    best_binding_site_score = bin_site_score
                    best_lin_score = lin_score
                    best_rmsd = rmsd
                    best_id = mutant_id
                line = sf.readline()
        return ScoreResult(best_score,
                           best_id + '.pdb',
                           binder_score=best_binder_score,
                           ligand_score=best_ligand_score,
                           binding_site_score=best_binding_site_score,
                           lin_score=best_lin_score,
                           rmsd_score=best_rmsd)

    def score(self, structure_file: str, descriptor: ComplexDescriptor, score_file_: str = None, log_file: str = None,
              delete_existing_results: bool = True) -> ScoreResult:
        """Invokes Rosetta relax in a shell, waits for the termination of the program and returns the result.
        :param structure_file: the path of a PDB file to relax
        :param descriptor: defines ligand and binder residues
        :param score_file_: path to a file where the output is stored
        :param log_file: Rosetta's command line output is directed into this file
        :param delete_existing_results: whether structures originating from previous Rosetta runs (probably based on
        entirely different mutation sets due to randomization!) shall be removed
        :return: the output of the relax, including the total score
        """
        if score_file_ is None:
            score_file_ = structure_file.lower().replace(".pdb", ".score")
        if log_file is None:
            log_file = structure_file.lower().replace(".pdb", ".log")
        output_path = self.mutant_dir
        with open(log_file, 'w') as lf:
            if delete_existing_results:
                for f in glob(structure_file.replace('.pdb', '_*.pdb')):
                    os.remove(f)
                sf = structure_file.replace('.pdb', '.score')
                if os.path.exists(sf):
                    os.remove(sf)
            executable = self.relax_executable
            if self.mode in ('score', 'fixbb'):
                executable = executable.replace('relax', self.mode)
            run_command(rosetta_path=self.rosetta_path, exe_file=executable, db_path=self.rosetta_db_path,
                        input_file=structure_file, output_dir=output_path, score_file=score_file_,
                        extra_rotamers=self.extra_rotamers, mode=self.mode, nstruct=self.nstruct,
                        script=self.relax_script, stdout=lf)
        if self.mode == 'score':
            shutil.copyfile(f"{structure_file}", f"{structure_file.replace('.pdb', '_0001.pdb')}")
        return self.read_score(score_file_, descriptor)

    def switch_to_nd(self, back: bool = False):
        """
        Alters mutant directory in negative design mode (and back)
        :param back: Reverse to normal mode.
        :return: None
        """
        if not back:
            with open("resfile.res", "w") as rf:
                rf.write("NATAA\nstart")
            self.mutant_dir = self.mutant_dir + "negative_design/"
            self.mode = 'fixbb'
            self.extra_rotamers = True
        else:
            if os.path.exists("./resfile.res"):
                os.remove("./resfile.res")
            self.mutant_dir = self.mutant_dir.split("negative_design/")[0]
            self.mode = 'relax'
