r"""ABC Index descriptor.

References
    * http://match.pmf.kg.ac.rs/electronic_versions/Match75/n1/match75n1_233-242.pdf

"""

import numpy as np

from ._base import Descriptor
from ._graph_matrix import DistanceMatrix

__all__ = ("ABCIndex", "ABCGGIndex")


class ABCIndexBase(Descriptor):
    __slots__ = ()

    @classmethod
    def preset(cls, version):
        yield cls()

    explicit_hydrogens = False

    def parameters(self):
        return ()

    def __str__(self):
        return self.__class__.__name__[:-5]

    rtype = float


class ABCIndex(ABCIndexBase):
    r"""atom-bond connectivity index descriptor.

    References:
        * :doi:`10.2298/FIL1204733D`

    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "atom-bond connectivity index"

    @staticmethod
    def _each_bond(bond):
        du = bond.GetBeginAtom().GetDegree()
        dv = bond.GetEndAtom().GetDegree()

        return np.sqrt((du + dv - 2.0) / (du * dv))

    def calculate(self):
        return float(sum(self._each_bond(bond) for bond in self.mol.GetBonds()))


class ABCGGIndex(ABCIndexBase):
    r"""Graovac-Ghorbani atom-bond connectivity index descriptor.

    References:
        * Furtula, B. Atom-bond connectivity index versus Graovac-Ghorbani analog. MATCH Commun. Math. Comput. Chem 75, 233-242 (2016).

    """  # noqa: E501

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "Graovac-Ghorbani atom-bond connectivity index"

    def dependencies(self):
        return {"D": DistanceMatrix(self.explicit_hydrogens)}

    @staticmethod
    def _each_bond(bond, D):
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()

        nu = np.sum(D[u, :] < D[v, :])
        nv = np.sum(D[v, :] < D[u, :])

        return np.sqrt((nu + nv - 2.0) / (nu * nv))

    def calculate(self, D):
        return float(sum(self._each_bond(bond, D) for bond in self.mol.GetBonds()))
