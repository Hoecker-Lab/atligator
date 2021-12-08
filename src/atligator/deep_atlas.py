"""This module contains facilities for deep learning atlases and predicting pmatrices from acomplexes based on a
convolutional neural network implemented by SciKit Learn.

:Deprecated: This module has been deprecated because machine learning did not show any improvement compared to the
traditional atlas solution, which implements in fact a linear classification with deterministic parametrization.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-03-14
"""

from pickle import dump, load
from typing import List, Dict, Set, Tuple

from Bio.PDB.vectors import Vector

from atligator.acomplex import LigandBinderComplex
from atligator.atlas import Atlas
from atligator.pdb_util import canonical_amino_acids
from atligator.prediction_matrix import BinderPredictionMatrix, BinderPredictionLine, BinderPredictionColumn
from atligator.structure import DesignedBinderResidue, DesignedLigandResidue, get_icoor_from_aresidue


def get_labels_from_resname(restype: str) -> List[float]:
    """Converts a residue type into a feature vector of length 20, where all positions contain the defined negativity
    value, 0.0, except for the index of the residue in the canonical list, which will contain the positivity value 1.0.
    :param restype: three letter code of the residue type to encode
    :returns: a vector of floats with length 20 that represents the feature vector that can be fed into a classifier
    """
    vec = [0.0] * len(canonical_amino_acids)
    try:
        index = canonical_amino_acids.index(restype)
        vec[index] = 1.0
    except ValueError:
        pass
    return vec


def get_labels_from_vector(vector: Vector, scale: float = 1.0) -> List[float]:
    """Scales a vector in order to better fit into an MLP classifier.
    :param vector: the vector to scale. Coordinates in Angstrom
    :param scale: the scalar value to multiply the vector with. If unset, no scaling is applied (factor 1.0)
    :return: a list with three floats that represent the scaled vector that can be fed to MLP
    """
    return (vector ** scale).get_array().tolist()


class DeepAtlas:
    """Instances of this class hold an MLP classifier incorporating the training data extracted from an atlas. The
    underlying classification problem is understood as follows:
    Ligand Residue Type x Ligand C_alpha x Ligand C_alpha->C_beta x Ligand C_alpha->C_O ==> Binder Residue Type
    The topology of the neural network is as follows:
    Input Layer: 29 neurons (20 for ligand residue type label, 3 each for binder C_alpha, orientation and secor coors)
    Hidden Layer: 25 neurons
    Output Layer 20 neurons (representing the possible values for the binder residue type label)"""

    def __init__(self, classifier=None, max_distance: float = 8.0, hidden_neurons: int = 25,
                 activation: str = 'tanh', learning_rate: float = 0.005, max_iter: int = 500, alpha: float = 0.01,
                 tolerance: float = 1e-6, verbose_training: bool = False):
        """Creates a new instance.
        :param classifier: if set, this will be the classifier of the deep atlas. Otherwise, a new MLP classifier is
        created based on experience parameters.
        :param max_distance: the maximum distance between C_alpha atoms in the underlying atlas (used for scaling
        the input vectors)
        :param hidden_neurons: size of the hidden layer of the underlying neural classificator. According to Jeff
        Heaton's rule of thumb, pick a layer size that is the mean of input and output neurons. Therefore, the default
        is 25.
        :param activation: the activation function to use for backpropagation. The default, tanh, works well with
        continuous inputs between -1.0 and 1.0, but also with discrete inputs between 0.0 and 1.0
        :param learning_rate: how fast (but also how biased) learning will happen
        :param max_iter: maximum number of iterations per training
        :param alpha: this factor trains the sharpness of the output
        :param tolerance: learning is stopped when training loss is below this value
        :param verbose_training: whether to pass over klearn's console output
        """
        if classifier is not None:
            self.classifier = classifier
        else:
            from sklearn.neural_network import MLPClassifier
            self.classifier = MLPClassifier(
                hidden_layer_sizes=(hidden_neurons,),
                activation=activation,
                learning_rate_init=learning_rate,
                max_iter=max_iter,
                alpha=alpha,
                tol=tolerance,
                verbose=verbose_training
            )
        self.max_distance = max_distance

    def train(self, atlas: Atlas) -> None:
        """Trains the classifier with input data originating from an atlas.
        :param atlas: its datapoints will be used as training samples.
        """
        data: List[List[float]] = []
        targets: List[List[float]] = []
        for ares in atlas.datapoints:
            sample = get_labels_from_resname(ares.ligand_restype)
            sample.extend(get_labels_from_vector(ares.calpha, 1.0 / self.max_distance))
            sample.extend(get_labels_from_vector(ares.orientation().normalized()))
            sample.extend(get_labels_from_vector(ares.secondary_orientation().normalized()))
            data.append(sample)
            targets.append(get_labels_from_resname(ares.binder_restype))
        self.classifier.fit(data, targets)

    def predict(self, binder_res: DesignedBinderResidue, ligand_res: DesignedLigandResidue,
                prob_threshold: float = 0.001) -> Dict[str, float]:
        """Using the underlying classificator, this method predicts the residue type of a of a given binder residue,
        relative to a given ligand residue, to whose interal coordinate system all binder vectors are transformed first.
        :param binder_res: the binder residue whose residue type to predict
        :param ligand_res: the context ligand residue.
        :param prob_threshold: classification results with probabilities below this threshold will be discarded
        :param: a dictionary that maps significantly probable residue types to the predicted probability
        """
        icoor = get_icoor_from_aresidue(ligand_res)
        sample = get_labels_from_resname(ligand_res.restype)
        sample.extend(get_labels_from_vector(icoor.external_to_internal(
            binder_res.calpha(), False), 1.0 / self.max_distance))
        sample.extend(get_labels_from_vector(icoor.external_to_internal(
            binder_res.orientation(), True).normalized()))
        sample.extend(get_labels_from_vector(icoor.external_to_internal(
            binder_res.secondary_orientation(), True).normalized()))
        prediction: List[float] = self.classifier.predict_proba([sample])[0]
        result: Dict[str, float] = {}
        for i, v in enumerate(prediction):
            restype: str = canonical_amino_acids[i]
            value = prediction[i]
            if value > prob_threshold:
                result[restype] = float(value)
        return result

    def save(self, filename) -> None:
        """Saves the classificator in binary format in the specified file.
        :param filename: where to serialize"""
        with open(filename, 'wb') as fo:
            dump(self, fo)


def load_deep_atlas(filename: str) -> DeepAtlas:
    """Creates a new atlas based on a serialized classificator instance.
    :param filename: from where to load the classificator
    :return: a corresponding deep atlas
    """
    with open(filename, 'rb') as fo:
        deep_atlas: DeepAtlas = load(fo)
    return deep_atlas


def generate_deep_prediction_matrix(acomplex: LigandBinderComplex, deep_atlas: DeepAtlas, distance_cutoff: float = 8.0,
                                    distance_exponent: float = 2.0) -> BinderPredictionMatrix:
    """Predicts possible residue types of the designed residues of a given ligand binder complex based on a deep atlas.
    To this end, every pair of designed ligand and binder residue, which lie within a given cutoff radius, is fed to
    the classificator. The local predictions of all ligand residues are superimposed and stored in the resulting
    binder prediction matrix
    :param acomplex: the ligand binder complex to create the matrix for
    :param deep_atlas: the deep atlas that is used for classification-based prediction
    :param distance_cutoff: maximum distance between ligand and binder C_alpha atoms
    :param distance_exponent: The relative distance (between 0.0 and 1.0) will be raised with this value before being
            multiplied to the score. Therefore, this exponent denotes the "steepness" of distance influence
    :return: a binder prediction matrix calculated based on the classification as explained
    """
    lines: List[BinderPredictionLine] = []
    for binder_res in acomplex.binder.residues:
        columns_map: Dict[str, Tuple[str, Set[int], float]] = {}
        for ligrespos, ligand_res in enumerate(acomplex.ligand.residues):
            distance: float = (ligand_res.calpha - binder_res.calpha).norm()
            if distance < distance_cutoff:
                prediction = deep_atlas.predict(binder_res, ligand_res)
                for pred_restype, pred_score in prediction.items():
                    score = float(pred_score) * (((distance_cutoff - distance) / distance_cutoff) ** distance_exponent)
                    if pred_restype not in columns_map:
                        columns_map[pred_restype] = (pred_restype, set(), 0.0)
                    col = columns_map[pred_restype]
                    columns_map[pred_restype] = (col[0], col[1].union({ligrespos}), col[2]+score)
        columns: List[BinderPredictionColumn] = []
        for col in columns_map.values():
            columns.append(BinderPredictionColumn(col[0], len(col[1]), col[2]))
        lines.append(BinderPredictionLine(binder_res.residue_id, binder_res.original_residue_type, columns))
    return BinderPredictionMatrix(lines)
