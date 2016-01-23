from ._base import Descriptor
from ._common import DistanceMatrix

import rdkit.Chem.GraphDescriptors as RDKit


class BertzCT(Descriptor):
    r"""Bertz CT descriptor(rdkit wrapper).

    :rtype: float
    """

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return 'BertzCT'

    def dependencies(self):
        return dict(D=DistanceMatrix(self.explicit_hydrogens))

    def calculate(self, mol, D):
        return RDKit.BertzCT(mol, dMat=D)
