"""The outcome of ATLIGATOR is the prediction of probable mutations on the binder site, increasing the binding
between binder and ligand. This module provides data structures and functions for mutation scoring. Scoring results
are kept in a so called binder prediction matrix.

:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-21
"""

from copy import deepcopy
from typing import List, Dict, Tuple, Set

from numpy import pi

from atligator.acomplex import LigandBinderComplex
from atligator.structure import AbstractResidue


class BinderPredictionColumn:
    """A column of a binder prediction line aggregates a residue type (indicating a mutation), a certainty value
    (number of atlas residues indicating the mutation, and the score value."""

    def __init__(self, residue_type: str, certainty: int, score: float):
        """Creates a new instance.
        :param residue_type: the candidate target type of a mutation
        :param certainty: the number of ligand residues from whose local atlases the mutation has been inferred
        :param score: relative probability of mutation success based on the overall atlas data
        """
        self.residue_type = residue_type
        self.certainty = certainty
        self.score = score

    def quality(self) -> float:
        """The quality of a column entry is defined as the product of score and certainty.
        :return: quality of this column
        """
        return self.score * float(self.certainty)

    def merge(self, other: 'BinderPredictionColumn') -> 'BinderPredictionColumn':
        """Merges this column with another column from another matrix by addding their scores and taking the maximum
        of their certainties.
        :param other: the column to merge with
        :return: a merged column
        """
        return BinderPredictionColumn(self.residue_type, max(self.certainty, other.certainty), self.score + other.score)


class BinderPredictionLine:
    """A line in a binder prediction aggregates the following data: the residue id of the involved binder residue,
    the original residue type, and a list of columns, each representing a suggested mutation."""

    def __init__(self, residue_id: int, original_restype: str, columns: List[BinderPredictionColumn]):
        """Creates a new instance.
        :param residue_id: the id of the binder residue candidate to mutation
        :param original_restype: the original residue type according to the scaffold
        :param columns: a list of columns, each predicting possible mutations for the context residue
        """
        self.residue_id = residue_id
        self.original_restype = original_restype
        self.columns = columns

    def merge(self, other: 'BinderPredictionLine') -> 'BinderPredictionLine':
        """Merges the columns of this matrix with corresponding columns of another, given matrix. As matching
        criteria for columns, we use the residue type.
        :param other: the line to merge with
        :return: a merged line containing merged columns
        """
        column_restypes = sorted(set([c.residue_type for c in self.columns + other.columns]))
        merged_columns: List[BinderPredictionColumn] = []
        for cr in column_restypes:
            if any([sc for sc in self.columns if sc.residue_type == cr]):
                sc = [sc for sc in self.columns if sc.residue_type == cr][0]
                if any([oc for oc in other.columns if oc.residue_type == cr]):
                    oc = [oc for oc in other.columns if oc.residue_type == cr][0]
                    merged_columns.append(sc.merge(oc))
                else:
                    merged_columns.append(deepcopy(sc))
            elif any([oc for oc in other.columns if oc.residue_type == cr]):
                oc = [oc for oc in other.columns if oc.residue_type == cr][0]
                merged_columns.append(deepcopy(oc))
        return BinderPredictionLine(self.residue_id, self.original_restype, merged_columns)


def read_binder_prediction_line(line: str) -> BinderPredictionLine:
    """Unparses a single line of a binder prediction matrix.
    :param line: the text serialization
    :return: a corresponding instance of BinderPredictionLine, including its columns
    """
    content = line.split()
    orig_restype = content[0]
    binder_res_id = int(content[1])
    columns: List[BinderPredictionColumn] = []
    for i in range(3, len(content), 3):
        residue_type = content[i]
        certainty = int(content[i+1])
        score = float(content[i+2])
        columns.append(BinderPredictionColumn(residue_type, certainty, score))
    return BinderPredictionLine(binder_res_id, orig_restype, columns)


class BinderPredictionMatrix:
    """A binder affinity matrix assigns to each designed residue of a ligand-binder complex a set of predicted
    scores for residue types."""

    def __init__(self, lines: List[BinderPredictionLine]):
        """Creates a new BinderPredictionMatrix from its lines.
        :param lines: The initial contents of the matrix.
        """
        self.lines = lines

    def __mul__(self, scalar: float) -> 'BinderPredictionMatrix':
        """Multiplies the scores of this matrix with a given scalar and returns the result
        :param scalar: the scalar value to multiply the score of each matrix cell with
        :return: a scaled prediction matrix
        """
        result_lines: List[BinderPredictionLine] = []
        for line in self.lines:
            result_columns: List[BinderPredictionColumn] = []
            for col in line.columns:
                result_columns.append(BinderPredictionColumn(col.residue_type, col.certainty, col.score * scalar))
            result_lines.append(BinderPredictionLine(line.residue_id, line.original_restype, result_columns))
        return BinderPredictionMatrix(result_lines)

    def get_max_sumscore(self, criterion: str = 'score') -> float:
        """Returns the max_sumscore, which denotes the inverse of the factor which is applied for normalization.
        :param criterion: the normalization criterion. Must be one of 'score', 'certainty', or 'quality'
        :return: the maximum sum of scores of all matrix lines
        """
        if criterion not in ['score', 'certainty', 'quality']:
            pass
        max_sumscore = 0.0
        for line in self.lines:
            sum_scores = 0.0
            for col in line.columns:
                if criterion == 'score':
                    sum_scores += col.score
                elif criterion == 'certainty':
                    sum_scores += col.certainty
                elif criterion == 'quality':
                    sum_scores += col.quality()
            max_sumscore = max(max_sumscore, sum_scores)
        if max_sumscore == 0.0:
            max_sumscore = 1.0
        return max_sumscore

    def get_total_sumscore(self, criterion: str = 'score') -> float:
        """Returns the sum of all score, certainty, or quality values contained in this matrix.
        :param criterion: the normalization criterion. Must be one of 'score', 'certainty', or 'quality'
        :return: the sum of values of all matrix lines
        """
        if criterion not in ['score', 'certainty', 'quality']:
            pass
        total_sumscore = 0.0
        for line in self.lines:
            for col in line.columns:
                if criterion == 'score':
                    total_sumscore += col.score
                elif criterion == 'certainty':
                    total_sumscore += col.certainty
                elif criterion == 'quality':
                    total_sumscore += col.quality()
        return total_sumscore

    def total_quality(self, ignore_identities: bool = False) -> float:
        """
        :return: the sum of quality values of all columns in all lines of this matrix
        """
        total_quality = 0.0
        for line in self.lines:
            for col in line.columns:
                if ignore_identities and line.original_restype == col.residue_type:
                    continue
                total_quality += col.quality()
        return total_quality

    def normalize(self, criterion: str = 'score', total: bool = False) -> 'BinderPredictionMatrix':
        """Rescales the matrix such that the maximum sum of scores, certainty, or quality contained in the columns of
        any line equals one. Keeps the context matrix untouched and returns the normalized variant.
        :param criterion: the normalization criterion. Must be one of 'score', 'certainty', 'quality'
        :param total: if True, the total score/certain/quality rather than the maximum is taken for scale
        :return: a normalized matrix
        """
        if criterion not in ['score', 'certainty', 'quality']:
            pass
        if total:
            return self * (1.0 / self.get_total_sumscore(criterion))
        else:
            return self * (1.0 / self.get_max_sumscore(criterion))

    def normalize_quality(self) -> 'BinderPredictionMatrix':
        """Shorthand for normalizing by quality.
        :return: a matrix normalized by quality
        """
        return self.normalize('quality')

    def sort_columns(self) -> 'BinderPredictionMatrix':
        """Sorts the columns of every line in the matrix by descending score. Returns the sorted matrix.
        :return: a column-sorted matrix
        """
        sorted_lines: List[BinderPredictionLine] = []
        for line in self.lines:
            sorted_columns = sorted(line.columns, key=(lambda col: col.quality()), reverse=True)
            sorted_lines.append(BinderPredictionLine(line.residue_id, line.original_restype, sorted_columns))
        return BinderPredictionMatrix(sorted_lines)

    def sort_lines(self) -> 'BinderPredictionMatrix':
        """Sorts the lines of this matrix by binder residue ID. Returns the sorted matrix.
        :return: a matrix sorted by lines
        """
        return BinderPredictionMatrix(sorted(self.lines, key=lambda l: int(l.residue_id)))

    def merge(self, other: 'BinderPredictionMatrix') -> 'BinderPredictionMatrix':
        """Merges the lines of this matrix with the lines of another given matrix and returns a merged matrix. The
        matching criterion is the residue id represented by each line.
        :param other: the other matrix to merge with
        :return: the merged matrix
        """
        line_resids = sorted(set([l.residue_id for l in self.lines + other.lines]))
        merged_lines: List[BinderPredictionLine] = []
        for lr in line_resids:
            if any([sl for sl in self.lines if sl.residue_id == lr]):
                sl = [sl for sl in self.lines if sl.residue_id == lr][0]
                if any([ol for ol in other.lines if ol.residue_id == lr]):
                    ol = [ol for ol in other.lines if ol.residue_id == lr][0]
                    merged_lines.append(sl.merge(ol))
                else:
                    merged_lines.append(deepcopy(sl))
            elif any([ol for ol in other.lines if ol.residue_id == lr]):
                ol = [ol for ol in other.lines if ol.residue_id == lr][0]
                merged_lines.append(deepcopy(ol))
        return BinderPredictionMatrix(merged_lines)

    def print_matrix(self, score_threshold: float = 0.001) -> None:
        """Prints out the columns of the matrix line by line
        :param score_threshold: score values below this threshold will not be printed
        """
        for line in self.lines:
            s = f"{line.original_restype} {line.residue_id} : "
            for col in line.columns:
                if col.score > score_threshold:
                    s += f"{col.residue_type} {col.certainty} {'{:10.6f}'.format(col.score)} "
            print(s)

    def save(self, filename, score_threshold: float = 0.001):
        """Saves this pmatrix into the specified file.
        :param filename: the desired serialization location
        :param score_threshold: columns indicating scores below this parameter are ignored for serialization
        """
        with open(filename, 'w') as fo:
            for line in self.lines:
                s = f"{line.original_restype} {line.residue_id} : "
                for col in line.columns:
                    if col.score > score_threshold:
                        s += f"{col.residue_type} {col.certainty} {'{:10.6f}'.format(col.score)} "
                fo.write(f"{s}\n")


def read_matrix(filename: str) -> BinderPredictionMatrix:
    """Converts a file containing a serialized binder prediction matrix into a runtime instance.
    :param filename: path of the file to open for deserialization
    :return: a parsed instance of BinderPredictionMatrix
    """
    lines: List[BinderPredictionLine] = []
    with open(filename, 'r') as fo:
        for line in fo:
            lines.append(read_binder_prediction_line(line))
    return BinderPredictionMatrix(lines)


def proximity_score(binder_res: AbstractResidue, aligned_res: AbstractResidue, score_threshold: float = 6.0,
                    distance_factor: float = 1.0, orient_factor: float = 6.0/pi, secor_factor: float = 6.0/pi,
                    score_exponent: float = 2.0) -> float:
    """Calculates the proximity score between a designed binder residue and an aligned residue. The higher this score
    is, the 'more proximate' the residues are. The score is in the interval [0; 1] and is indirectly proportional to a
    weighted sum of the distance, the angle between c_alpha->c_beta, and the angle between c_alpha->c_o.
    :param binder_res: the designed binder residue to compare
    :param aligned_res: the aligned residue to compare
    :param score_threshold: if the weighted sum of deviation values is above this value, the score is constantly 0.0
    :param distance_factor: the share of C_alpha-C_alpha distance in the scoring function [1 / Angstrom]
    :param orient_factor: the share of the angle between orientation vectors [1 / rad]
    :param secor_factor: the share of the secondary orientation angle in the weighted sum [1 / rad]
    :param score_exponent: before returning, the calculated score is raised with this value
    :return: a score in the interval between [0; 1] that describes the result of structural comparison
    """
    distance: float = (aligned_res.calpha() - binder_res.calpha()).norm()
    orient_angle: float = binder_res.orientation().angle(aligned_res.orientation())
    sec_orient_angle: float = binder_res.secondary_orientation().angle(aligned_res.secondary_orientation())
    wsum = distance_factor * distance + orient_factor * orient_angle + secor_factor * sec_orient_angle
    if wsum > score_threshold:
        return 0.0
    else:
        score = ((score_threshold - wsum) / score_threshold) ** score_exponent
        return score


def generate_prediction_matrix(acomplex: LigandBinderComplex, **score_args) -> BinderPredictionMatrix:
    """Evaluates the given ligand-binder complex and calculates binder prediction matrix by applying the scoring
    function for every pair of binder residue and aligned residue.
    :param acomplex: the ligand-binder complex to use as input for the calculations
    :param score_args: additional arguments that will be passed to the scoring function. See proximity_score
    :return: a non-normalized and unsorted binder prediction matrix generated
    """
    lines: List[BinderPredictionLine] = []
    n_aligned_res = len(acomplex.aligned_residues)
    for binres in acomplex.binder.residues:
        res_id = binres.residue_id
        columns_map: Dict[str, Tuple[str, Set[int], float]] = {}
        for ares in acomplex.aligned_residues:
            score = proximity_score(binres, ares, **score_args) / float(ares.size_factor) / float(n_aligned_res)
            if score > 0.0:
                if ares.binder_restype not in columns_map:
                    columns_map[ares.binder_restype] = (ares.binder_restype, set(), 0.0)
                col: Tuple[str, Set[int], float] = columns_map[ares.binder_restype]
                columns_map[ares.binder_restype] = (col[0], col[1].union({ares.ligand_res_position}), col[2]+score)
        columns: List[BinderPredictionColumn] = []
        for col in columns_map.values():
            columns.append(BinderPredictionColumn(col[0], len(col[1]), col[2]))
        lines.append(BinderPredictionLine(res_id, binres.original_residue_type, columns))
    return BinderPredictionMatrix(lines)
