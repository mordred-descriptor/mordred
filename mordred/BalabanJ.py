import rdkit.Chem.GraphDescriptors as RDKit

from ._base import Descriptor
from ._graph_matrix import DistanceMatrix


__all__ = ('BalabanJ',)


class BalabanJ(Descriptor):
    r"""Balaban's J index descriptor(rdkit wrapper)."""

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def as_key(self):
        return self.__class__, ()

    def __str__(self):
        return 'BalabanJ'

    def dependencies(self):
        return {'D': DistanceMatrix(self.explicit_hydrogens)}

    def calculate(self, D):
        return float(RDKit.BalabanJ(self.mol, dMat=D))

    rtype = float
