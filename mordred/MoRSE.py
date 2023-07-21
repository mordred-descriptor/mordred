from __future__ import division

from itertools import chain

import numpy as np

from ._graph_matrix import DistanceMatrix3D
from ._atomic_property import AtomicProperty
from ._base.descriptor import Descriptor

__all__ = ("MoRSE",)


class MoRSE(Descriptor):
    since = "1.0.0"
    __slots__ = ("_prop", "_distance")

    require_3D = True

    def description(self):
        return "3D-MoRSE{} (distance = {})".format(
            ""
            if self._prop is None
            else " weighted by {}".format(self._prop.get_long()),
            self._distance,
        )

    @classmethod
    def preset(cls, version):
        return chain(
            (cls(None, i) for i in range(1, 33)),
            (cls("m", i) for i in range(1, 33)),
            (cls("v", i) for i in range(1, 33)),
            (cls("se", i) for i in range(1, 33)),
            (cls("p", i) for i in range(1, 33)),
        )

    def __str__(self):
        p = "" if self._prop is None else self._prop.as_argument
        return "Mor{:02d}{}".format(self._distance, p)

    def parameters(self):
        p = None if self._prop is None else self._prop.as_argument
        return p, self._distance

    def __init__(self, prop=None, distance=2):
        if prop is None:
            self._prop = None
        else:
            self._prop = AtomicProperty(self.explicit_hydrogens, prop)

        self._distance = distance

    def dependencies(self):
        d = {"D": DistanceMatrix3D(self.explicit_hydrogens)}

        if self._prop is not None:
            d["A"] = self._prop

        return d

    def calculate(self, D, A=None):
        if D.shape[0] <= 1:
            self.fail(ValueError("require 2 or more atoms"))

        N = D.shape[0]

        if A is None:
            A = np.ones(N)
        else:
            A = A / self._prop.carbon

        A = A.reshape(1, -1)

        if self._distance == 1:
            n = np.ones((N, N), dtype="float")

        else:
            with self.rethrow_zerodiv():
                sr = (self._distance - 1) * D
                np.fill_diagonal(sr, 1)
                n = np.sin(sr) / sr

        np.fill_diagonal(n, 0)

        return float(0.5 * np.ravel(A.dot(n).dot(A.T))[0])

    rtype = float
