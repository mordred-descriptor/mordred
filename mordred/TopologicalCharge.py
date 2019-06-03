from itertools import chain

import numpy as np
from six import integer_types

from ._base import Descriptor
from ._graph_matrix import DistanceMatrix, AdjacencyMatrix

__all__ = ("TopologicalCharge",)


class ChargeTermMatrix(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False

    def parameters(self):
        return ()

    def dependencies(self):
        return {
            "A": AdjacencyMatrix(self.explicit_hydrogens),
            "D": DistanceMatrix(self.explicit_hydrogens),
        }

    def calculate(self, A, D):
        D2 = D.copy()
        D2[D2 != 0] **= -2
        np.fill_diagonal(D2, 0)

        M = A.dot(D2)
        return M - M.T


class TopologicalCharge(Descriptor):
    r"""topological charge descriptor.

    :type type: str
    :param type:
        * "raw": sum of order-distance atom pairs coefficient
        * "mean": mean of order-distance atom pairs coefficient
        * "global": sum of mean-topoCharge over 0 to order

    :type order: int
    :param order: int

    References
        * :doi:`10.1021/ci00019a008`

    """

    since = "1.0.0"
    __slots__ = ("_type", "_order")

    explicit_hydrogens = False

    tc_types = ("global", "mean", "raw")

    def description(self):
        return "{}-ordered {} topological charge".format(self._order, self._type)

    @classmethod
    def preset(cls, version):
        return chain(
            (cls(t, o) for t in ("raw", "mean") for o in range(1, 11)),
            [cls("global", 10)],
        )

    def __str__(self):
        if self._type == "global":
            return "JGT{}".format(self._order)
        elif self._type == "mean":
            return "JGI{}".format(self._order)
        else:
            return "GGI{}".format(self._order)

    def parameters(self):
        return self._type, self._order

    def __init__(self, type="global", order=10):
        assert type in self.tc_types
        assert type == "global" or isinstance(order, integer_types)

        self._type = type
        self._order = order

    def dependencies(self):
        return {"CT": ChargeTermMatrix(), "D": DistanceMatrix(self.explicit_hydrogens)}

    def calculate(self, CT, D):
        D = D * np.tri(*D.shape)
        D[D == 0] = np.inf

        f = D <= self._order if self._type == "global" else D == self._order

        CT = CT[f]

        if self._type == "raw":
            return np.abs(CT).sum()

        # create frequency vector
        Df = D[f]
        C = Df.copy()
        for i in np.unique(Df):
            C[Df == i] = len(Df[Df == i])

        return np.abs(CT / C).sum()

    rtype = float
