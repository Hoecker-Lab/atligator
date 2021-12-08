from argparse import FileType, _ArgumentGroup, ArgumentParser
from typing import Dict

# Every dict needs a type. Type bool will only be taken for identification, not used in argparse.
# Every int, float and str type dict need a default value. str and FileType('r') defaults can be "", i.e. no default.
# Every str type dict needs choices. empty string ("") means there are no given choices.
# 'alter' key is used for alternative but not mandatory choices and is only considered if choices is "".
atlas_dict = {"group": "None", "type": "FileType('r')"}
scaffold_dict = {"group": "None", "type": "FileType('r')"}
descriptor_dict = {"group": "None", "type": "FileType('r')"}
g_dict = {"group": 'opt None', "l_name": 'graphical_input', "type": 'bool', "action": "store_true"}
a_dict = {"group": 'opt None', "l_name": 'atpar', "type": 'str', "default": "", "choices": ""}
lb_dict = {"group": 'lc', "l_name": 'max_siblings', "type": 'int', "default": 4}
lc_dict = {"group": 'lc', "l_name": 'n_clones', "type": 'int', "default": 1}
le_dict = {"group": 'lc', "l_name": 'n_epochs', "type": 'int', "default": 10}
lg_dict = {"group": 'lc', "l_name": 'max_generations', "type": 'int', "default": 10}
lg1_dict = {"group": 'lc', "l_name": 'gen_first_epoch', "type": 'int', "default": -1}
lk_dict = {"group": 'lc', "l_name": 'n_workers', "type": 'int', "default": -1, "limit": (-1, 128)}
ln_dict = {"group": 'lc', "l_name": 'n_top', "type": 'int', "default": 20}
lr_dict = {"group": 'lc', "l_name": 'n_recomb', "type": 'int', "default": 8}
fs_dict = {"group": 'fit', "l_name": 'scorer', "type": 'str', "default": 'rosetta',
           "choices": ["rosetta", "osprey", "simulate"]}
mc_dict = {"group": 'mut', "l_name": 'clustering_factor', "type": 'float', "default": 1.0}
md_dict = {"group": 'mut', "l_name": 'deep_factor', "type": 'float', "default": 0.0}
mf_dict = {"group": 'mut', "l_name": 'clash_factor', "type": 'float', "default": 1.0}
mr_dict = {"group": 'mut', "l_name": 'random_mutations', "type": 'bool', "action": "store_true"}
mrc_dict = {"group": 'mut', "l_name": 'random_mutation_chance', "type": 'int', "default": 20, "limit": (0, 100)}
mu_dict = {"group": 'mut', "l_name": 'coupling', "type": 'bool', "action": "store_true"}
sd_dict = {"group": 'sc', "l_name": 'distance_factor', "type": 'float', "default": 1.0}
se_dict = {"group": 'sc', "l_name": 'score_exponent', "type": 'float', "default": 2.0}
so_dict = {"group": 'sc', "l_name": 'orient_factor', "type": 'float', "default": 2.0}
ss_dict = {"group": 'sc', "l_name": 'secor_factor', "type": 'float', "default": 2.0}
st_dict = {"group": 'sc', "l_name": 'score_threshold', "type": 'float', "default": 8.0}
cc_dict = {"group": 'cs', "l_name": 'clash_merge_base', "type": 'float', "default": 1.5}
cd_dict = {"group": 'cs', "l_name": 'delta_factor', "type": 'float', "default": 4.0}
cn_dict = {"group": 'cs', "l_name": 'n_delta', "type": 'float', "default": 5.0}
cs_dict = {"group": 'cs', "l_name": 'std_dev_weight', "type": 'float', "default": 1.5}
cq_dict = {"group": 'cs', "l_name": 'n_quality_k', "type": 'float', "default": 5.0}
da_dict = {"group": 'da', "l_name": 'deep_alpha', "type": 'float', "default": 0.01}
dc_dict = {"group": 'da', "l_name": 'deep_activation', "type": 'str', "default": 'tanh',
           "choices": ['identity', 'logistic', 'tanh', 'relu']}
dd_dict = {"group": 'da', "l_name": 'deep_distance_cutoff', "type": 'float', "default": 8.0}
dh_dict = {"group": 'da', "l_name": 'deep_hidden_neurons', "type": 'int', "default": 25}
dl_dict = {"group": 'da', "l_name": 'deep_learning_rate', "type": 'float', "default": 0.005}
dm_dict = {"group": 'da', "l_name": 'deep_max_iter', "type": 'int', "default": 500}
dt_dict = {"group": 'da', "l_name": 'deep_tolerance', "type": 'float', "default": 1e-6}
dv_dict = {"group": 'da', "l_name": 'deep_verbose', "type": 'bool', "action": "store_true"}
rb_dict = {"group": 'rr', "l_name": 'rosetta_binding_site_factor', "type": 'float', "default": -1.0, "limit": (-1, 100)}
rd_dict = {"group": 'rr', "l_name": 'rosetta_db_path', "type": 'str', "default": "$ROSETTA_DB", "choices": ""}
re_dict = {"group": 'rr', "l_name": 'rosetta_extra_rotamers', "type": 'bool', "action": "store_true"}
rl_dict = {"group": 'rr', "l_name": 'rosetta_ligand_factor', "type": 'float', "default": 10.0}
rm_dict = {"group": 'rr', "l_name": 'rosetta_mode', "type": 'str', "default": 'relax',
           "choices": ['relax', 'score']}
rn_dict = {"group": 'rr', "l_name": 'rosetta_nstruct', "type": 'int', "default": 1}
ro_dict = {"group": 'rr', "l_name": 'rosetta_rmsd_offset', "type": 'float', "default": 3.0}
rp_dict = {"group": 'rr', "l_name": 'rosetta_path', "type": 'str', "default": "$ROSETTA", "choices": ""}
rq_dict = {"group": 'rr', "l_name": 'rosetta_pen_factor', "type": 'float', "default": 10.0}
rr_dict = {"group": 'rr', "l_name": 'rosetta_rmsd_slope', "type": 'float', "default": 2.0}
rs_dict = {"group": 'rr', "l_name": 'rosetta_relax_script', "type": "FileType('r')", "default": ""}
rv_dict = {"group": 'rr', "l_name": 'rosetta_negative_design', "type": 'bool', "action": "store_true"}
rx_dict = {"group": 'rr', "l_name": 'rosetta_relax_executable', "type": 'str',
           "default": 'relax.static.linuxgccrelease', "choices": ""}
ry_dict = {"group": 'rr', "l_name": 'rosetta_linearity_factor', "type": 'float', "default": 5.0}
pa_dict = {"group": 'os', "l_name": 'osprey_astar_mode', "type": 'str', "default": 'mplp',
           "choices": ['traditional', 'mplp']}
pe_dict = {"group": 'os', "l_name": 'osprey_use_epic', "type": 'bool', "action": "store_true"}
pf_dict = {"group": 'os', "l_name": 'osprey_finder', "type": 'str', "default": 'gmec', "choices": ['gmec', 'deegmec']}
pl_dict = {"group": 'os', "l_name": 'osprey_use_lute', "type": 'bool', "action": "store_true"}
ps_dict = {"group": 'os', "l_name": 'osprey_sidechain_flexibility', "type": 'str', "default": 'continuous',
           "choices": ['library', 'discrete', 'continuous']}
od_dict = {"group": 'out', "l_name": 'mutator_dumpfile', "type": 'str', "choices": "", "default": "",
           "alter": 'genmut.dump'}
ol_dict = {"group": 'out', "l_name": 'log_file', "type": "FileType('w')", "default": "genmut.log"}
om_dict = {"group": 'out', "l_name": 'mutant_dir', "type": 'str', "default": './mutants/', "choices": "", "alter": './'}
oo_dict = {"group": 'out', "l_name": 'output_file', "type": "FileType('w')", "default": 'genmut.score'}
op_dict = {"group": 'out', "l_name": 'mutant_pattern', "type": 'str', "default": 'mutant-$.pdb', "choices": ""}

arg: Dict[str, Dict] = {"atlas": atlas_dict,
                        "scaffold": scaffold_dict,
                        "descriptor": descriptor_dict,
                        "G": g_dict,
                        "A": a_dict,
                        "Lb": lb_dict,
                        "Lc": lc_dict,
                        "Le": le_dict,
                        "Lg": lg_dict,
                        "Lg1": lg1_dict,
                        "Lk": lk_dict,
                        "Ln": ln_dict,
                        "Lr": lr_dict,
                        "Fs": fs_dict,
                        "Mc": mc_dict,
                        "Md": md_dict,
                        "Mf": mf_dict,
                        "Mr": mr_dict,
                        "Mrc": mrc_dict,
                        "Mu": mu_dict,
                        "Sd": sd_dict,
                        "Se": se_dict,
                        "So": so_dict,
                        "Ss": ss_dict,
                        "St": st_dict,
                        "Cc": cc_dict,
                        "Cd": cd_dict,
                        "Cn": cn_dict,
                        "Cs": cs_dict,
                        "Cq": cq_dict,
                        "Da": da_dict,
                        "Dc": dc_dict,
                        "Dd": dd_dict,
                        "Dh": dh_dict,
                        "Dl": dl_dict,
                        "Dm": dm_dict,
                        "Dt": dt_dict,
                        "Dv": dv_dict,
                        "Rb": rb_dict,
                        "Rd": rd_dict,
                        "Re": re_dict,
                        "Rl": rl_dict,
                        "Rm": rm_dict,
                        "Rn": rn_dict,
                        "Ro": ro_dict,
                        "Rp": rp_dict,
                        "Rq": rq_dict,
                        "Rr": rr_dict,
                        "Rs": rs_dict,
                        "Rv": rv_dict,
                        "Rx": rx_dict,
                        "Ry": ry_dict,
                        "Pa": pa_dict,
                        "Pe": pe_dict,
                        "Pf": pf_dict,
                        "Pl": pl_dict,
                        "Ps": ps_dict,
                        "Od": od_dict,
                        "Ol": ol_dict,
                        "Om": om_dict,
                        "Oo": oo_dict,
                        "Op": op_dict}

help_txts: Dict[str, str] = \
    {"atlas": "the .atlas file to take as basis for genetic mutation.",
     "scaffold": "the PDB scaffold to optimize",
     "descriptor": "the descriptor file defining the ligand residue types as well as mutable host residues",
     "G": "Opens graphical interface for input arguments.",
     "A": "Load values from '.atpar' file. Only enabled in combination with graphical interface.",
     "Lb": "the maximum number of individuals with identical mutation sets that may coexist per epoch. Should be a "
           "multiple of n_clones",
     "Lc": "this number of clones is created for every individual generated. Clones are then processedin parallel, "
           "each is carrying a unique individual ID. Original is also counted.",
     "Le": "the number of epochs of genetic optimization",
     "Lg": "the number of recombining generations per epoch",
     "Lg1": "the number of recombining generations per epoch in first epoch.",
     "Lk": "Number of workers for parallel processing recombinations. -1 for number of CPU cores",
     "Ln": "the number of top elements to compare after each iteration. If the ranking obtained this way does not "
           "change, the script is aborted. Top mutations are also printed out regularly.",
     "Lr": "the number best mutants to recombine with base mutations in each generation",
     "Fs": "the scorer to use for calculating the fitness and the optimized structure of individuals",
     "Mc": "probability of distance clustering based mutation selection (default selection strategy)",
     "Md": "probability to involve a deep learning strategy for mutation selection",
     "Mf": "probability for clash factoring of the selection of mutations in each selection round",
     "Mr": "whether to apply random mutations.",
     "Mrc": "The chance of random mutations (over ones based on pmatrix). The chance is given in percent and only "
            "taken into account, if random mutations are enabled.",
     "Mu": "whether to apply analogous mutations to coupled residues automatically",
     "Sd": "weight factor for distance coefficient of the scoring function. [1/Angstrom]",
     "Se": "the outcome of the scoring function is postprocessed by raising by this value.",
     "So": "weight factor for the primary orientation coefficient (C_alpha->C_beta) of the scoring function. "
           "[1/radians]",
     "Ss": "weight factor for the secondary orientation coefficient (C_alpha->C_O) of the scoring function. "
           "[1/radians]",
     "St": "threshold coefficient of the scoring function for residue alignment during the pmatrix calculation "
           "in every epoch.",
     "Cc": "Base of clash punishment for merging with normalized quality: n_quality*BASE^-clash.",
     "Cd": "Correction factor for merging clash/quality with delta within calculation of mutation "
           "ranking: 'n_quality*base^-clash + n_delta/FACTOR'.",
     "Cn": "Normalization constant for delta. Equals to delta value of 50%% in normalized state.",
     "Cs": "Base in weighing the standard deviation: BASE^-std_dev * clash.",
     "Cq": "Normalization constant for quality. Equals to quality value of 50%% in normalized state.",
     "Da": "a coefficient that controls how sharp classification output will become (lower = sharper)",
     "Dc": "the activation function to use for the backpropagation algorithm while learning.",
     "Dd": "the minimum length of a chain to be treated as ligand",
     "Dh": "the size of the hidden layer in the underlying neural classificator.",
     "Dl": "learning rate of the neural classificator.",
     "Dm": "maximum number of iterations to process for each training data set",
     "Dt": "a coefficient that controls how small the training effect must be to terminate tranining",
     "Dv": "verbose training output",
     "Rb": "The share of binding site energy in the calculation of total energy by Rosetta relax. If -1, it is set"
           "to the ligand factor value.",
     "Rd": "Path or environment variable pointing to the Rosetta database folder",
     "Re": "Whether Rosetta relax uses extra rotamers (Rosetta flags -ex1 -ex2 set)",
     "Rl": "The share of ligand energy in the calculation of total energy by Rosetta relax",
     "Rm": "Rosetta's execution mode. Use 'simulate' for testing purposes only. 'score' will not relax the "
           "structure, but just score it, which is much faster yet unrealistic.",
     "Rn": "number of starting points for every Rosetta relax run. Hint: consider parameter n_clones instead, "
           "which behaves almost equivalently yet allows for parallel Rosetta runs.",
     "Ro": "An offset value where rmsd penalty will start",
     "Rp": "Path or environment variable pointing to the folder containing Rosetta executables",
     "Rq": "weight factor of the RMSD penalty",
     "Rr": "The slope of the influence of ligand rmsd after the offset",
     "Rs": "Rosetta relax script to use. If unset, use default relax:fast protocol. Ignored if rosetta_mode "
           "is 'simulate'.",
     "Rv": "Use negative design for best individuals to test specificity",
     "Rx": "Name of the relax executable in rosetta_path",
     "Ry": "Impact of peptide linearity on the calculation of total energy by Rosetta relax",
     "Pa": "A* implementation",
     "Pe": "whether to use the EPIC algorithm for continuous flexibility modeling",
     "Pf": "Finder method",
     "Pl": "whether to use the LUTE algorithm for continuous flexibility modeling",
     "Ps": "Resolution for sidechain optimization",
     "Od": "if this file exists, the genetic mutator will be reconstructed from it. Furthermore, after every "
           "generation, the current progress is dumped into this file.",
     "Ol": "if set, the log output is forwarded to this file, otherwise to stdout",
     "Om": "directory where mutant pdb and log files are stored. Will be created if not existant.",
     "Oo": "output file. Linewise output in format <total-score> <epoch> <generation> <mutant-id> <mutations>, "
           "ordered by ascending total-score.",
     "Op": "file name pattern for mutated structures. '$' will be replaced by the line number of the corresponding "
           "mutation in the mutations file."}


def add_argument(group: ArgumentParser or _ArgumentGroup, name: str, optional: bool = True):
    if optional:
        arg_name = "-" + name
        l_name = "--" + arg[name]['l_name']
        if arg[name]['type'] == 'bool':
            group.add_argument(arg_name, l_name, action=arg[name]['action'], help=help_txts[name])
        elif arg[name]['type'] == 'int':
            group.add_argument(arg_name, l_name, default=arg[name]['default'], type=int,
                               help=help_txts[name])
        elif arg[name]['type'] == 'float':
            group.add_argument(arg_name, l_name, default=arg[name]['default'], type=float,
                               help=help_txts[name])
        elif arg[name]['type'] == 'str':
            if arg[name]['choices'] != "":
                if arg[name]['default'] != "":
                    group.add_argument(arg_name, l_name, default=arg[name]['default'], type=str,
                                       choices=arg[name]['choices'], help=help_txts[name])
                else:
                    group.add_argument(arg_name, l_name, type=str,
                                       choices=arg[name]['choices'], help=help_txts[name])
            else:
                if arg[name]['default'] != "":
                    group.add_argument(arg_name, l_name, default=arg[name]['default'], type=str,
                                       help=help_txts[name])
                else:
                    group.add_argument(arg_name, l_name, type=str, help=help_txts[name])
        elif arg[name]['type'] == "FileType('r')":
            if arg[name]['default'] != "":
                group.add_argument(arg_name, l_name, default=arg[name]['default'], type=FileType('r'),
                                   help=help_txts[name])
            else:
                group.add_argument(arg_name, l_name, type=FileType('r'), help=help_txts[name])
        elif arg[name]['type'] == "FileType('w')":
            if arg[name]['default'] != "":
                group.add_argument(arg_name, l_name, default=arg[name]['default'], type=FileType('w'),
                                   help=help_txts[name])
            else:
                group.add_argument(arg_name, l_name, type=FileType('w'), help=help_txts[name])

    else:
        if arg[name]['type'] == 'bool':
            group.add_argument(name, action=arg[name]['action'], help=help_txts[name])
        elif arg[name]['type'] == 'int':
            group.add_argument(name, type=int, help=help_txts[name])
        elif arg[name]['type'] == 'float':
            group.add_argument(name, type=float, help=help_txts[name])
        elif arg[name]['type'] == 'str':
            if arg[name]['choices'] != "":
                group.add_argument(name, type=str, choices=arg[name]['choices'], help=help_txts[name])
            else:
                group.add_argument(name, type=str, help=help_txts[name])
        elif arg[name]['type'] == "FileType('r')":
            group.add_argument(name, type=FileType('r'), help=help_txts[name])
        elif arg[name]['type'] == "FileType('w')":
            group.add_argument(name, type=FileType('w'), help=help_txts[name])
    return group
