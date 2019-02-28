import numpy as np

from ._base import Descriptor
from ._graph_matrix import DistanceMatrix3D
from ._atomic_property import AtomicProperty, get_properties

__all__ = ("RDF",)


def _distances():
    for i in range(2, 32):
        yield 0.5 * i


class RDF(Descriptor):
    __slots__ = ("_prop", "_R", "_beta")

    require_3D = True

    since = "1.1.0"

    def description(self):
        return "RDF descriptor(beta={}) - {}. weighted by {}".format(
            self._beta,
            self._R,
            self._prop.get_long(),
        )

    @classmethod
    def preset(cls, version):
        for p in get_properties(version):
            for r in _distances():
                yield cls(r, p)

    def __str__(self):
        name = "RBF{:03}{}".format(int(self._R * 10), self._prop.as_argument)

        if self._beta != 100:
            name += "/{}".format(self._beta)

        return name

    def parameters(self):
        return (self._R, self._prop, self._beta)

    def __init__(self, R=1.0, prop="x", beta=100):
        self._prop = AtomicProperty(self.explicit_hydrogens, prop)
        self._R = R
        self._beta = beta

    def dependencies(self):
        return {
            "D": DistanceMatrix3D(True),
            "P": self._prop,
        }

    def calculate(self, P, D):
        Pr = P / self._prop.carbon
        W = Pr[:, np.newaxis].dot(Pr[np.newaxis, :])

        e = np.exp(-self._beta * (self._R - D) ** 2)
        e *= W
        np.fill_diagonal(e, 0)
        return np.sum(e) / 2

    rtype = float
