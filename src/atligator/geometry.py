"""Data structures and operations for handling geometric representations relative to different internal or external
coordinate systems.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2018-03-05
"""

from typing import List

from Bio.PDB.vectors import Vector
from numpy.linalg import inv

# typedefs
Matrix = List[List[float]]


class InternalCoordinates:
    """Instances of this class represent an internal coordinate system, which is defined by means of a translation
    vector and a rotation matrix."""

    def __init__(self, translation_vector: Vector, rotation_matrix: Matrix):
        """creates an instance of an internal coordinate system. It is defined by a translation vector and a rotation
        matrix.
        :param translation_vector vector between external and internal origin
        :param rotation_matrix matrix describing rotation from external to internal coordinates. This matrix should
        be normalized.
        """
        self.translation_vector = translation_vector
        self.rotation_matrix = rotation_matrix

    def external_to_internal(self, external_vector: Vector, direction: bool) -> Vector:
        """transforms a given external vector into a corresponding vector relative to the internal coordinate system.
        :param external_vector: vector in external coordinates to convert
        :param direction: if false, both translation and rotation will be applied; otherwise, only rotation
        :return: the vector in internal coordinates
        """
        if direction:
            internal_vector = external_vector
        else:
            internal_vector = external_vector - self.translation_vector
        internal_vector = internal_vector.left_multiply(self.rotation_matrix)
        return internal_vector

    def internal_to_external(self, internal_vector: Vector, direction: bool) -> Vector:
        """transforms a given vector relative to this internal coordinate system into a corresponding external
        (absolute) vector.
        :param internal_vector: vector to convert in internal coordinates
        :param direction: if true, both rotation and translation will be applied; otherwise, only rotation
        :return: the vector in external coordinates
        """
        external_vector = internal_vector.left_multiply(inv(self.rotation_matrix))
        if not direction:
            external_vector = external_vector + self.translation_vector
        return external_vector
