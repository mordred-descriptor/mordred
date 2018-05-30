import rdkit.Chem.GraphDescriptors as RDKit

from ._base import Descriptor
from ._graph_matrix import DistanceMatrix

__all__ = ("BertzCT",)


class BertzCT(Descriptor):
    r"""Bertz CT descriptor(rdkit wrapper)."""

    since = "1.0.0"
    __slots__ = ()
    explicit_hydrogens = False

    def description(self):
        return "Bertz CT"

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return self.__class__.__name__

    def parameters(self):
        return ()

    def dependencies(self):
        return {"D": DistanceMatrix(self.explicit_hydrogens)}

    def calculate(self, D):
        return float(RDKit.BertzCT(self.mol, dMat=D))

    rtype = float
