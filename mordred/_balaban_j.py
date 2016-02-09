import rdkit.Chem.GraphDescriptors as RDKit

from ._base import Descriptor
from ._graph_matrix import DistanceMatrix


class BalabanJ(Descriptor):
    r"""Balaban's J index descriptor(rdkit wrapper)."""

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def __reduce_ex__(self, version):
        return self.__class__, ()

    def __str__(self):
        return 'BalabanJ'

    def dependencies(self):
        return {'D': DistanceMatrix(self.explicit_hydrogens)}

    def calculate(self, mol, D):
        return float(RDKit.BalabanJ(mol, dMat=D))

    rtype = float
