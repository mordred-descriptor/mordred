import numpy as np

from ._base import Descriptor
from ._util import atoms_to_numpy
from ._graph_matrix import AdjacencyMatrix, DistanceMatrix3D

__all__ = ("GravitationalIndex",)


class GravitationalIndex(Descriptor):
    since = "1.0.0"
    __slots__ = "_heavy", "_pair"

    def description(self):
        return "{}{}gravitational index".format(
            "heavy atom " if self._heavy else "", "pair " if self._pair else ""
        )

    @classmethod
    def preset(cls, version):
        return (cls(h, p) for p in [False, True] for h in [True, False])

    require_3D = True

    @property
    def explicit_hydrogens(self):
        return not self._heavy

    def __str__(self):
        return "GRAV{}{}".format("" if self._heavy else "H", "p" if self._pair else "")

    def parameters(self):
        return self._heavy, self._pair

    def __init__(self, heavy=True, pair=False):
        self._heavy = heavy
        self._pair = pair

    def dependencies(self):
        d = {"D": DistanceMatrix3D(self.explicit_hydrogens)}

        if self._pair:
            d["A"] = AdjacencyMatrix(self.explicit_hydrogens)

        return d

    def calculate(self, D, A=1.0):
        w = atoms_to_numpy(lambda a: a.GetMass(), self.mol)

        w = w[:, np.newaxis] * w
        np.fill_diagonal(w, 0)

        D = D.copy()
        np.fill_diagonal(D, 1)

        with self.rethrow_zerodiv():
            return 0.5 * np.sum(w * A / D ** 2)

    rtype = float
