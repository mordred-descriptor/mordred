from ._base import Descriptor
from ._common import DistanceMatrix

import rdkit.Chem.GraphDescriptors as RDKit


class BalabanJ(Descriptor):
    r"""Balaban's J index descriptor(rdkit wrapper).

    :rtype: float
    """

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return 'BalabanJ'

    def dependencies(self):
        return dict(D=DistanceMatrix(self.explicit_hydrogens))

    def calculate(self, mol, D):
        return RDKit.BalabanJ(mol, dMat=D)
