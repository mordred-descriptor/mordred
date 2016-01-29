import math

import numpy as np

from ._base import Descriptor
from ._common import DistanceMatrix


class ABCIndexBase(Descriptor):
    @classmethod
    def preset(cls):
        yield cls()

    explicit_hydrogens = False

    def __reduce_ex__(self, version):
        return self.__class__, ()

    def __str__(self):
        return self.__class__.__name__[:-5]

    rtype = float


class ABCIndex(ABCIndexBase):
    r"""atom-bond connectivity indez descriptor.
    """

    def calculate(self, mol):
        s = 0.0

        for bond in mol.GetBonds():
            du = bond.GetBeginAtom().GetDegree()
            dv = bond.GetEndAtom().GetDegree()

            s += math.sqrt(float(du + dv - 2) / (du * dv))

        return s


class ABCGGIndex(ABCIndexBase):
    r"""Graovac-Ghorbani atom-bond connectivity index descriptor.
    """

    def dependencies(self):
        return {'D': DistanceMatrix(self.explicit_hydrogens)}

    def calculate(self, mol, D):
        s = 0.0

        for bond in mol.GetBonds():
            u = bond.GetBeginAtomIdx()
            v = bond.GetEndAtomIdx()

            nu = np.sum(D[u, :] < D[v, :])
            nv = np.sum(D[v, :] < D[u, :])

            s += np.sqrt(float(nu + nv - 2) / (nu * nv))

        return s
