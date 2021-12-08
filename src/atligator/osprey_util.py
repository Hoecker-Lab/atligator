"""This module provides an implementation of the BaseScorer interface that utilizes the Osprey protein design suite.
Unfortunately, the Osprey python interface seems to be strictly limited to python2. Therefore, we use another
indirection here: our python3 module creates a python2 script, which is executed and starts the Osprey Java app.

For installing Osprey we used a local installation ('python2 get-pip.py --user' from 'https://bootstrap.pypa.io/get-pip.py') of pip (pip2! NOT inside conda env!).
After downloading (https://github.com/donaldlab/OSPREY_refactor/releases) and unzipping the osprey package you have
to run following commands:
    'pip2 install osprey --no-index --user --find-link=wheelhouse'
    'pip2 install -U numpy --user'


THIS recommended version will NOT work: (wheel is deprecated and installation has to be done locally)
    'pip2 install osprey --no-index --use-wheel --find-link=wheelhouse'


:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-04-19
"""

from subprocess import run

from atligator.acomplex import read_complex_descriptor
from atligator.base_scorer import BaseScorer, ScoreResult
from atligator.genetic_mutator import ComplexDescriptor
from atligator.pdb_util import replace_restypes


class OspreyScorer(BaseScorer):
    """Implements the BaseScorer interface upon Osprey."""

    def __init__(self, sidechain_flexibility: str, astar_mode: str, finder: str, use_epic: bool, use_lute: bool,
                 python2_executable: str = "python2", descriptor: ComplexDescriptor = None, cores: int = None):
        """Creates a new instance.
        :param sidechain_flexibility: Resolution for sidechain optimization. Must be one of ['library', 'discrete',
        'continuous']
        :param astar_mode: A* implementation to use. Must be one of ['traditional', 'mplp]
        :param finder: Finder method to use for optimization. Must be one of ['gmec', 'deegmec']
        :param use_epic: whether to use the EPIC algorithm for continuous flexibility modeling
        :param use_lute: whether to use the LUTE algorithm
        :param python2_executable: the path to python2 to pass to the shell
        :param descriptor: path to descriptor file to read out mutable binder and ligand residues
        :param cores: number of cores used for energy calculation
        """
        if sidechain_flexibility not in ['library', 'discrete', 'continuous', 'wildtype']:  # 'wildtype' for rigid
            raise ValueError()
        self.sidechain_flexibility = sidechain_flexibility
        if astar_mode not in ['traditional', 'mplp']:
            raise ValueError()
        self.astar_mode = astar_mode
        if finder not in ['gmec', 'deegmec', 'lute']:  # 'lute' for new lute_train method and testing
            raise ValueError()
        self.finder = finder
        self.use_epic = use_epic
        self.use_lute = use_lute
        self.python2_executable = python2_executable
        if descriptor is None:
            self.descriptor = read_complex_descriptor("./descriptor")
        else:
            self.descriptor = descriptor
        self.cores = 1 if cores is None else cores

    def score(self, structure: str, descriptor: ComplexDescriptor) -> ScoreResult:
        """
        Optimizes and scores the structure according to the Osprey parameters set in the constructor.
        :param structure: structure to optimize and score
        :param descriptor: defines ligand and binder residues
        :return: optimized structure and score obtained
        """
        replace_restypes(structure, 'HIS', 'HID')
        optimized_structure = structure.replace('.pdb', '_opt.pdb')
        res_bin = []
        res_lig = []
        for i in self.descriptor.binder_info:
            res_bin.append(i[0] + i[1])
        for i in self.descriptor.ligand_info:
            res_lig.append(i[0] + i[1])

        print(res_lig, res_bin)
        # Last few lines are needed to read out which residues to set flexible (mutatable ones from descriptor)
        # To toggle between these and all residues change commenting in the lines:
        # f"for res in strand.mol.residues:\n" + \
        # f"#for resname in {res_bin}:\n" + \
        # f"    resname = res.pDBResNumber\n"
        py2script = f"import osprey\n" + \
                    f"osprey.start()\n" + \
                    f"strand = osprey.Strand('{structure}')\n" + \
                    f"#for res in strand.mol.residues:\n" + \
                    f"for resname in {res_bin+res_lig}:\n" + \
                    f"#    resname = res.pDBResNumber\n"
        # .setLibraryRotamers(osprey.WILD_TYPE) sets to any residue type specified in the brackets. In this case
        # a library of rotamers of the wild type residue type (aa) is used. empty brackets would give same result.
        # .addWildTypeRotamers() adds the wild type rotamer. So with 'wildtype' only the original rotamer (or rotamers
        # if there are several conformations) is used.
        if self.sidechain_flexibility == 'library':
            py2script += f"    strand.flexibility[resname].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()\n"
        elif self.sidechain_flexibility == 'discrete':
            py2script += f"    strand.flexibility[resname].setLibraryRotamers(osprey.WILD_TYPE).setDiscrete()\n"
        elif self.sidechain_flexibility == 'wildtype':
            py2script += f"    strand.flexibility[resname].addWildTypeRotamers()\n"
        else:
            py2script += f"    strand.flexibility[resname].setLibraryRotamers(osprey.WILD_TYPE).setContinuous()\n"
        # Not sure if ff_params.solvScale = 0. is necessary. Found in one of x examples.
        #            f"bbflex = osprey.DEEPerStrandFlex(strand, '1CC8.pert', {res_lig}, '{structure}')\n" \
        #            f"conf_space = osprey.ConfSpace([[strand, bbflex]])\n" + \
        py2script += f"conf_space = osprey.ConfSpace(strand)\n" + \
                     f"ff_params = osprey.ForcefieldParams()\n" + \
                     f"ff_params.solvScale = 0.\n" + \
                     f"ene_calc = osprey.EnergyCalculator(conf_space, ff_params, " \
                     f"parallelism=osprey.Parallelism(cpuCores={self.cores}))\n" + \
                     f"conf_calc = osprey.ConfEnergyCalculator(conf_space, ene_calc)\n" + \
                     f"ene_mat = osprey.EnergyMatrix(conf_calc)\n"
        if self.astar_mode == 'mplp':
            py2script += f"astar_alg = osprey.AStarMPLP(ene_mat, conf_space)\n"
        else:
            py2script += f"astar_alg = osprey.AStarTraditional(ene_mat, conf_space)\n"
        if self.finder == 'gmec':
            py2script += f"finder = osprey.GMECFinder(astar_alg, conf_calc)\n" + \
                         f"solution = finder.find()\n"
        else:
            py2script += f"finder = osprey.DEEGMECFinder(ene_mat, conf_space, ene_calc, conf_calc, 'deegmec', " + \
                         f"use_epic={self.use_epic}, use_lute={self.use_lute})\n" + \
                         f"solution = finder.calcGMEC()\n"
        py2script += f"score = solution.getScore()\n" + \
                     f"optimized_molecule = conf_space.makeMolecule(solution.getAssignments())\n" + \
                     f"osprey.writePdb(optimized_molecule, '{optimized_structure}')\n" + \
                     f"print score\n" + \
                     f"print '{optimized_structure}'\n"
        # Step to use new lute_gmecFinder method of osprey. Not implemented in current version (3.0b1), but found in
        # examples/documentation and also found in git repository. If new osprey version gets released, uncomment this.
        # 71+72+75+78+79+81+82+83+88+113+114+117+120+121+123+124+125+130+155+156+159+162+163+165+166+167+172+197+198+201+204+205+207
        """
        if self.finder == 'lute':
            py2script = f"import osprey\n" + \
                        f"osprey.start()\n" + \
                        f"strand = osprey.Strand('{structure}')\n" \
                        f"for res in strand.mol.residues:\n" + \
                        f"    resname = res.pDBResNumber\n" \
                        f"    strand.flexibility[resname].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()\n" \
                        f"confSpace = osprey.ConfSpace(strand)\n" \
                        f"pmat = osprey.DEE_read(confSpace, 'LUTE.pmat.dat')\n" \
                        f"model = osprey.LUTE_read('LUTE.dat')\n" \
                        f"gmecFinder = osprey.LUTE_GMECFinder(confSpace, model, pmat, printIntermediateConfs=True)\n" \
                        f"gmec = gmecFinder.find(1.0)\n" \
                        f"score = solution.getScore()\n" \
                        f"optimized_molecule = conf_space.makeMolecule(solution.getAssignments())\n" + \
                        f"osprey.writePdb(optimized_molecule, '{optimized_structure}')\n" + \
                        f"print score\n" + \
                        f"print '{optimized_structure}'\n"
        """
        py2file = structure.replace(".pdb", "-osprey.py")
        with open(py2file, 'w') as fo:
            fo.write(py2script)
        py2out = structure.replace(".pdb", "-osprey.out")
        with open(py2out, 'w') as pfo:
            run(f"{self.python2_executable} {py2file}", shell=True, check=True, stdout=pfo)
        with open(py2out, 'r') as pfo:
            lines = list(pfo.readlines())
            score = float(lines[-2])
            optimized_structure = lines[-1].strip()
        replace_restypes(structure, 'HID', 'HIS')
        replace_restypes(optimized_structure, 'HID', 'HIS')
        return ScoreResult(score, optimized_structure)
