"""This module provides an interface for different scoring functions that can be used, e.g., with the genetic mutator.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-04-19
"""

import random
from shutil import copyfile

from atligator.acomplex import ComplexDescriptor


class ScoreResult:
    """Instances of this class represent the result of a scoring done by a subclass of BaseScorer."""

    def __init__(self, score: float, optimized_structure: str, **kwargs):
        """Creates a new instance.
        :param score: the score calculated
        :param optimized_structure: the conformation that represents the solution with the best score
        :param kwargs: additional properties to assign to this score result
        """
        self.score = score
        self.optimized_structure = optimized_structure
        self.properties = kwargs

    def __getattr__(self, item):
        return self.properties[item]

    def __str__(self):
        return str(self.score) + ' [' + ', '.join([str(k) + ': ' + str(v) for k, v in self.properties.items()]) + ']'

    def __repr__(self):
        return str(self)


class BaseScorer:
    """This abstract base class provides an interface for different scorers. Concrete scorers must implement the
    score function."""

    def score(self, structure: str, descriptor: ComplexDescriptor) -> ScoreResult:
        """To be implemented by subclasses. Evaluates the structure and returns the resulting score and optimized
        structure
        :param structure: the structure to optimize
        :param descriptor: defines ligand and binder residues of the structure
        :return: the corresponding score result
        """
        raise NotImplementedError


class RandomScorer(BaseScorer):
    """This scorer returns a random value between -1000.0 and 0.0. Use it for testing purposes only."""

    def score(self, structure: str, descriptor: ComplexDescriptor) -> ScoreResult:
        relaxed_file = structure.replace(".pdb", "_0001.pdb")
        copyfile(structure, relaxed_file)
        score = random.uniform(-1000.0, 0.0)
        return ScoreResult(score, relaxed_file)
