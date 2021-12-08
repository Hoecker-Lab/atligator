"""

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-15
"""

import math
from typing import Dict, List, Tuple

from atligator.io_util import HiddenPrints, Colors
from atligator.mutation import Mutation
from atligator.prediction_matrix import BinderPredictionMatrix, BinderPredictionLine, BinderPredictionColumn

"""
Hydrophobicity of a residue.
Derived from: J. Kyte, R. F. Doolittle: A simple method for displaying the hydropathic character of a protein. 
Journal of Molecular Biology. Band 157, Nr. 1, 1982, S. 105–132, PMID 7108955. 
Recalculated from -10 (Arg: -4.5) to 10 (Ile: 4.5). 
"""
hp_index = {"ALA": 4.00, "ARG": -10.00, "ASN": -7.78, "ASP": -7.78, "CYS": 5.56,
            "GLN": -7.78, "GLU": -7.78, "GLY": -0.89, "HIS": -7.11, "ILE": 10.00,
            "LEU": 8.44, "LYS": -8.67, "MET": 4.22, "PHE": 6.22, "PRO": -3.56,
            "SER": -1.78, "THR": -1.56, "TRP": -2.00, "TYR": -2.89, "VAL": 9.33}


charge_index = {"ALA": 0.00, "ARG": 10.00, "ASN": 0.00, "ASP": -10.00, "CYS": 0.00,
                "GLN": 0.00, "GLU": -10.00, "GLY": 0.00, "HIS": 2.00, "ILE": 0.00,
                "LEU": 0.00, "LYS": 10.00, "MET": 0.00, "PHE": 0.00, "PRO": 0.00,
                "SER": 0.00, "THR": 0.00, "TRP": 0.00, "TYR": 0.00, "VAL": 0.00}


"""
Relative size of a residue
Derived from van der Waals radius. Recalculated from -10 (Gly: 48) to 10 (Trp: 163). 
"""
size_index = {"ALA": -6.70, "ARG": 7.39, "ASN": -1.65, "ASP": -2.52, "CYS": 3.39,
              "GLN": 1.48, "GLU": 0.61, "GLY": -10.00, "HIS": 2.17, "ILE": 3.22,
              "LEU": 3.22, "LYS": 5.13, "MET": 3.22, "PHE": 5.13, "PRO": -2.70,
              "SER": -5.65, "THR": -2.17, "TRP": 10.00, "TYR": 6.17, "VAL": -0.09}


def generate_clash_pmatrix(pmatrix: BinderPredictionMatrix,
                           parameters: Tuple[float, float, float, float, float]) -> BinderPredictionMatrix:
    max_scoresum = pmatrix.get_max_sumscore()
    pmatrix = pmatrix.normalize()
    with HiddenPrints():
        similarity_matrix = generate_similarity_scores(pmatrix, parameters, max_scoresum)
    return similarity_matrix.generate_derived_pmatrix()


class FeatureSet:
    """
    A set which combines one value for each feature (hydrophobicity, size and charge).
    """

    def __init__(self, hydrophobicity: float, size: float, charge: float, h_dev: float, s_dev: float, c_dev: float):
        self.hydrophobicity = hydrophobicity
        self.size = size
        self.charge = charge
        self.h_dev = h_dev
        self.s_dev = s_dev
        self.c_dev = c_dev

    def __str__(self) -> str:
        """converts a feature set into its string representation.
        :return: string serialization
        """
        hp = "{0:.2f}".format(self.hydrophobicity)
        h_dev = "{0:.2f}".format(self.h_dev)
        size = "{0:.2f}".format(self.size)
        s_dev = "{0:.2f}".format(self.s_dev)
        charge = "{0:.2f}".format(self.charge)
        c_dev = "{0:.2f}".format(self.c_dev)
        pm = " +-"
        string = "hp: " + hp + pm + h_dev + "   size: " + size + pm + s_dev + "   charge: " + charge + pm + c_dev
        return string

    def __repr__(self):
        return str(self)


class SimilarityMatrix:
    """
    A Matrix combining residue id (position) with all mutations.
    These mutations are combined by three letter code for amino acid (residue id), residue name, all mutations
    of this position with clash scores to average feature set and original residue (=delta) (generated
    by calculate_clash) and total quality
    for each position on the binding interface.
    """
    def __init__(self, positions: List,
                 max_scoresum: float,
                 base: float = 1.5,
                 delta_factor: float = 4.0,
                 n_delta_k: float = 5.0,
                 n_quality_k: float = 5.0):
        self.positions = positions
        self.n_delta_k = n_delta_k
        self.n_quality_k = n_quality_k
        self.base = base
        self.max_score = max_scoresum
        self.delta_factor = delta_factor

    def __iter__(self):
        for i in range(len(self.positions)):
            yield self.positions[i]

    def generate_derived_pmatrix(self) -> BinderPredictionMatrix:
        lines = []
        for position in self.positions:
            columns = []
            res_id = position[0]
            res_type = position[1]
            total_quality = position[3]
            for mutation in position[2]:
                mut_res_type = mutation[0]
                clash = mutation[1]
                delta = mutation[2]
                rating = self.rate(clash, total_quality, delta)
                columns.append(BinderPredictionColumn(mut_res_type, 1, rating))
            lines.append(BinderPredictionLine(res_id, res_type, columns))
        return BinderPredictionMatrix(lines)

    def rate(self, clash: float, total_quality: float, delta: float):
        n_delta = delta / (delta + self.n_delta_k)
        tq = total_quality
        clash_n_quality = self.base ** (-clash) * tq * self.max_score / (tq * self.max_score + self.n_quality_k)
        rating = n_delta/self.delta_factor + clash_n_quality
        return rating


def generate_similarity_scores(pmatrix: BinderPredictionMatrix, parameters: Tuple[float, float, float, float, float],
                               max_scoresum: float = 1.0) -> SimilarityMatrix:
    """
    Generates three mutations suggestions for each position in a binder prediction matrix.
    Takes P matrix as an imput where all mutations are read (similar to mutation_prediction.find_best_mutations)
    For each position an average value for all features (hydrophobicity, size and charge)
    gets calculated and compared with the values of all amino acids.
    The best three matches get returned
    :param pmatrix: a binder prediction matrix
    :param parameters: Parameters for weighing: clash_merge_base, delta_factor, n_delta, n_quality_k, std_dev_weight
    :param max_scoresum: The maximum sum of score (height of pmatrix bar) for one position
    :return: A SimilarityMatrix instance including clash scores
    """
    mutations = find_all_mutations(pmatrix)
    hydrophobicity_matrix = calculate_hydrophobicity(mutations)
    size_matrix = calculate_size(mutations)
    charge_matrix = calculate_charge(mutations)
    feature_matrix = []
    for i in range(len(size_matrix)):
        h = hydrophobicity_matrix[i][2]
        s = size_matrix[i][2]
        c = charge_matrix[i][2]
        h_dev = hydrophobicity_matrix[i][3]
        s_dev = size_matrix[i][3]
        c_dev = charge_matrix[i][3]
        features = FeatureSet(h, s, c, h_dev, s_dev, c_dev)
        residue_id = size_matrix[i][0]
        residue_name = size_matrix[i][1]
        total_quality = size_matrix[i][4]
        feature_matrix.append((residue_id, residue_name, features, total_quality))
    similarity_matrix = calculate_clash_matrix(feature_matrix, max_scoresum, parameters)
    return similarity_matrix


def find_all_mutations(pmatrix: BinderPredictionMatrix) -> List[List[Tuple[int, Mutation]]]:
    """
    Reads BinderPredictionMatrix object and returns all mutations with their quality.
    In contrast to find_best_mutations all mutations are considered - even exchanging to
    the same residue (= no mutation basically).
    These mutations are not sorted by their impact which is not necessary since their processed
    all the same way in feature calculation.
    :param pmatrix: Binder prediction matrix with mutations including score and certainty
    :return: List of all mutations which can be found in the pmatrix.
    """
    mutations = []
    for line in pmatrix.lines:
        position = []
        for col in line.columns:
            quality = col.score * int(col.certainty)
            mut = Mutation(line.residue_id, line.original_restype, col.residue_type)
            position.append((quality, mut))
        mutations.append(position)
    return mutations


def calculate_clash_matrix(feature_matrix: List[Tuple[int, str, FeatureSet, float]], max_scoresum: float,
                           par: Tuple[float, float, float, float, float]) -> SimilarityMatrix:
    """
    Calculates a matrix with all clash scores for all mutations in the p matrix,
    but just keeps the three mutations (at one position) with the lowest clash score.
    :param feature_matrix: Combines FeatureSets for all positions including resid, resname and total_quality
    :param max_scoresum: The maximum sum of score (height of pmatrix bar) for one residue.
    :param par: Parameters for weighing: clash_merge_base, delta_factor, n_delta, n_quality_k, std_dev_weight
    :return: A SimilarityMatrix object including clash scores
    """
    mutation_matrix = []
    for position in feature_matrix:
        res_id = position[0]
        res_name = position[1]
        total_quality = position[3]
        feature_set = position[2]
        lowest_clash = calculate_clash(feature_set, par[4], res_name)
        mutation_matrix.append((res_id, res_name, lowest_clash, total_quality))
    similarity_matrix = SimilarityMatrix(mutation_matrix, max_scoresum, par[0], par[1], par[2], par[3])
    return similarity_matrix


def calculate_clash_dep(features: FeatureSet, base: float = 1.5) -> Tuple[List, List, List]:
    """
    DEPRECATED
    Calculates a total clash for all amino acids regarding the current average feature values.
    Returns best three mutations.
    """
    best = [None, 30]
    second = [None, 30]
    third = [None, 30]
    for residue in hp_index:
        hp_clash = base ** - features.h_dev * abs(hp_index[residue] - features.hydrophobicity)
        size_clash = base ** - features.s_dev * abs(size_index[residue] - features.size)
        charge_clash = base ** - features.c_dev * abs(charge_index[residue] - features.charge)
        total_clash = hp_clash + size_clash + charge_clash
        if total_clash < best[1]:
            third = second
            second = best
            best = [residue, total_clash]
        elif total_clash < second[1]:
            third = second
            second = [residue, total_clash]
        elif total_clash < third[1]:
            third = [residue, total_clash]
    return best, second, third


def calculate_clash(features: FeatureSet, base: float, original_res: str = None) -> List[Tuple]:
    """
    Calculates a total clash for all amino acids regarding the current average feature values.
    The standard deviation allows to introduce a weight on different features:
        -> If standard deviation is high: Feature accuracy might be less important
        -> If standard deviation is low: "        "        might be very important
    base value for weighing can be adjusted (see parameters). If base == 1, no weighing is applied.

    If an original_res is defined, an additional clash is calculated:
    The clash between the current mutation residue and the original residue.
    This can serve as an indicator for a change in the features.
    ATTENTION: This calculation is also weighted, if base is not set to 1.

    Returns all legal (= stated in binder prediction matrix) mutations sorted by increasing clash score.
    """
    mutations = [('', 30.0, 0.0)]
    for residue in hp_index:
        hp_clash: float = base ** - features.h_dev * abs(hp_index[residue] - features.hydrophobicity)
        size_clash: float = base ** - features.s_dev * abs(size_index[residue] - features.size)
        charge_clash: float = base ** - features.c_dev * abs(charge_index[residue] - features.charge)
        total_clash: float = hp_clash + size_clash + charge_clash
        for mut in range(len(mutations)):
            if total_clash < mutations[mut][1]:
                clash_to_original = 0.0
                if original_res is not None:
                    hp_clash = base ** - features.h_dev * abs(hp_index[residue] - hp_index[original_res])
                    size_clash = base ** - features.s_dev * abs(size_index[residue] - size_index[original_res])
                    charge_clash = base ** - features.c_dev * abs(charge_index[residue] - charge_index[original_res])
                    clash_to_original = hp_clash + size_clash + charge_clash
                mutations.insert(mut, (residue, total_clash, clash_to_original))
                break
    mutations.pop()
    return mutations


def calculate_feature(index: Dict[str, float], mutations: List[List[Tuple[float, Mutation]]]) \
        -> List[Tuple[int, str, float, float, float]]:
    """
    Calculates the average score for each position.
    The type of score (feature: size, charge, etc.) is defined by the index which is passed.
    :param index: The index defining the feature type.
    :param mutations:
    :return:
    """
    feature_matrix = []
    counter = -1
    for position in mutations:
        counter += 1
        sum_ = 0
        total_quality = 0
        for mutation in position:
            feature_value = index[mutation[1].mutated_residue_type]
            quality = mutation[0]
            sum_ += quality * feature_value
            total_quality += quality
        if total_quality == 0:
            print("position", counter, "skipped due to missing data input.")
            continue
        average = sum_ / total_quality
        dev_sum = 0.0
        for mutation in position:
            feature_value = index[mutation[1].mutated_residue_type]
            quality = mutation[0]
            dev_sum += (feature_value - average)**2 * quality
        variance = dev_sum / total_quality
        std_dev = math.sqrt(variance)
        res_id = position[0][1].residue_id
        res_name = position[0][1].original_residue_type
        feature_matrix.append((res_id, res_name, average, std_dev, total_quality))
    return feature_matrix


def calculate_hydrophobicity(mutations):
    return calculate_feature(hp_index, mutations)


def calculate_size(mutations):
    return calculate_feature(size_index, mutations)


def calculate_charge(mutations):
    return calculate_feature(charge_index, mutations)


def print_output(similarity_matrix: SimilarityMatrix, args, pmatrix_max_sumscore, colors: Colors) -> None:
    """
    Prints similarity scores for each position.
    :param similarity_matrix: A matrix with residue id, residue name, mutations sorted by increasing clash
    and total quality for each position on the binding interface.
    :param args: Arguments given by argparser
    :param pmatrix_max_sumscore: The maximum sum of score (height of pmatrix bar) for one position
    :param colors: Colors instance for coloring
    :return: None
    """
    # Normalisation constant: k is the value of delta where it reaches 50 % in normalized state
    n_delta_k = args.n_delta_k
    n_quality_k = args.n_quality_k
    base = args.base
    max_ = pmatrix_max_sumscore
    print("\nMaximum score sum is:", max_)
    for position in similarity_matrix:
        # Intro: "ASN 45 - mutation:"
        if len(str(position[0])) == 3:
            print(position[1], position[0], "-           mutation", end=":   ")
        else:
            print(position[1], position[0], " -           mutation", end=":   ")
        # Header with residue type - also defines how many are displayed (counter):
        counter = 0
        for i in range(len(position[2])):
            if base ** (-position[2][i][1]) * position[3] * max_ / (position[3] * max_ + n_quality_k) < 0.05:
                break
            if position[2][i][0] == position[1]:
                print(colors.RED + position[2][i][0], " (o)  ", end=colors.WHITE)
            else:
                print(position[2][i][0], end="      ")
            counter += 1
            if counter == args.output_size_ex:
                break

        # clash scores - based on similarity of pmatrix data and static amino acid features - of best mutations
        # within range (atm controlled by normalized quality threshold)
        print("\n clash-score:  ┬─>to ori(del):   ", end="")

        if args.norm_del:
            for i in range(counter):
                delta = position[2][i][2]
                print("┼d>" + colors.DIM, "{0:.1f}".format(delta), " ", end=colors.WHITE)
            print(colors.WHITE + colors.BOLD + "\n clash-score:  │ └>norm_delta:   ", end="")
            for i in range(counter):
                delta = position[2][i][2]
                n_delta = delta / (delta + n_delta_k)
                print("└d>", "{0:.1f}".format(n_delta), " ", end="")
        else:
            colors.pc('bold')
            for i in range(counter):
                delta = position[2][i][2]
                print("┴d>", "{0:.1f}".format(delta), " ", end="")

        print(colors.WHITE + "\n               └─>to average :   ", end="")
        for i in range(counter):
            print("┬a>", "{0:.1f}".format(position[2][i][1]), " ", end="")

        # quality of mutation in pmatrix. Defined by score of clustering and certainty.
        # clash_quality is a mixed value by clash score and quality. defined by ' base^(-clash) * quality '
        if args.quality:
            print("\n  quality:  ", "{0:.2f}".format(position[3]), " - clash-qu:   ", end="")
            for i in range(counter):
                clash_quality = base ** (-position[2][i][1]) * position[3]
                print("├q>" + colors.DIM, "{0:.1f}".format(clash_quality), end="  " + colors.WHITE)

        # normalized quality. Takes into account that big differences at high qualities might not be worth punishing.
        # the basis is michaelis menten kinetics-like formula: 'value = (max_[=100%] * quality)/(quality + 50%-value)'
        # This leads to close to linear rating at small quality values and tiny differences at high quality values.
        norm_quality = 100 * position[3] * max_ / (position[3] * max_ + n_quality_k)
        norm_quality_str = "{0:.2f}".format(norm_quality) if norm_quality >= 10 else " " + "{0:.2f}".format(
            norm_quality)
        colors.pc('white')
        colors.pc('bold')
        print("\n n-quality: ", norm_quality_str, "% - clash-qu: ", sep="", end="  ")
        for i in range(counter):
            clash_quality = base ** (-position[2][i][1]) * norm_quality / 100
            print("└n>", "{0:.1f}".format(clash_quality), end="  ")
        print("\n" + colors.WHITE)


def rate_mutations(similarity_matrix: SimilarityMatrix, args, pmatrix_max_sumscore) -> \
        List[Tuple[str, str, str, float]]:
    """
    Sorts mutations in similarity matrix and returns list with list from best to worst.
    This list includes residue id (position), original residue type, mutated residue type and a rating.
    The rating is based on the formula: rating = norm_delta / delta_factor + clashed_quality.
    delta (the clash of mutation residue to original) is normalized at high values. This can be adjusted by
    n_delta_k (see argparse help).
    clashed_quality is calculated by base^-(clash to average) * norm_quality. base and quality normalization
    factor can be adjusted with argparse.
    :param similarity_matrix: A matrix with residue id, residue name, mutations sorted by increasing clash
    and total quality for each position on the binding interface.
    :param args: Arguments given by argparser
    :param pmatrix_max_sumscore: The maximum sum of score (height of pmatrix bar) for one position
    :return: None
    """
    sort = [("position", "original", "mutation", 0.0)]
    base = args.base
    # Normalisation constant: k is the value of delta where it reaches 50 % in normalized state
    n_delta_k = args.n_delta_k
    n_quality_k = args.n_quality_k
    # Correction factor for (n_quality*base^-clash + n_delta/FACTOR)
    delta_factor = args.delta_factor
    max_ = pmatrix_max_sumscore

    for position in similarity_matrix:
        # Header with residue type:
        for i in range(len(position[2])):
            if position[1] == position[2][i][0]:
                continue
            delta = position[2][i][2]
            n_delta = delta / (delta + n_delta_k)
            clash_n_quality = base ** (-position[2][i][1]) * position[3] * max_ / (position[3] * max_ + n_quality_k)
            for mut in range(len(sort)):
                # print(type(n_delta + clash_n_quality), type(sort[mut][2]))
                rating = n_delta / delta_factor + clash_n_quality
                if rating > sort[mut][3]:
                    sort.insert(mut, (position[0], position[1], position[2][i][0], rating))
                    break
    sort.pop()
    return sort


def show_mutations(similarity_matrix: SimilarityMatrix, args, pmatrix_max_sumscore, colors: Colors) -> None:
    """
    Showing the best 10 (or '-os') mutations in a human readable ranking.
    With the color flas '-c' all the doubled ones are dimmed after appearing lateron.
    :param similarity_matrix: A matrix with residue id, residue name, mutations sorted by increasing clash
    and total quality for each position on the binding interface.
    :param args: Arguments given by argparser
    :param pmatrix_max_sumscore: The maximum sum of score (height of pmatrix bar) for one position
    :param colors: Colors instance for coloring
    :return: None
    """
    sort = rate_mutations(similarity_matrix, args, pmatrix_max_sumscore)
    counter = 0
    print("rating original to mutation  (n_quality*base^-clash + n_delta/factor)")
    ids_out = []
    for i in sort:
        if i[0] in ids_out:
            colors.pc('dim')
        else:
            colors.pc('white')
            ids_out.append(i[0])
        if counter == args.output_size:
            break
        if len(str(counter + 1)) == 2:
            print(counter + 1, ": ", sep="", end="")
        else:
            print(counter + 1, ":  ", sep="", end="")
        print("  ", i[1], "", end="")
        dyn_space = ""
        for j in range(3 - len(str(i[0]))):
            dyn_space += " "
        print(i[0], dyn_space, sep="", end=" ")
        print(" to ", i[2], ":      {0:.2f}".format(i[3]), sep="")
        counter += 1


def safe_output(similarity_matrix: SimilarityMatrix, args, pmatrix_max_sumscore) -> None:
    """
    Saving a Mutation Set with the best 10 mutations (or changed by '-os' flag).
    :param similarity_matrix: A matrix with residue id, residue name, mutations sorted by increasing clash
    and total quality for each position on the binding interface.
    :param args: Arguments given by argparser
    :param pmatrix_max_sumscore: The maximum sum of score (height of pmatrix bar) for one position
    :return: None
    """
    sort = rate_mutations(similarity_matrix, args, pmatrix_max_sumscore)
    counter = 0

    for i in sort:
        if counter == args.output_size:
            break
        print(i[1], i[0], i[2], sep="", end=" ")
        counter += 1
