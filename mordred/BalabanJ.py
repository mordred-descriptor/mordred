import rdkit.Chem.GraphDescriptors as RDKit

from ._base import Descriptor
from ._graph_matrix import DistanceMatrix

__all__ = ("BalabanJ",)


class BalabanJ(Descriptor):
    r"""Balaban's J index descriptor(rdkit wrapper)."""

    __slots__ = ()

    def description(self):
        return "Balaban's J index"

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def parameters(self):
        return ()

    def __str__(self):
        return self.__class__.__name__

    def dependencies(self):
        return {"D": DistanceMatrix(self.explicit_hydrogens)}

    def calculate(self, D):
        return float(RDKit.BalabanJ(self.mol, dMat=D))

    rtype = float
