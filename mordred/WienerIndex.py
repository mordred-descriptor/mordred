from ._base import Descriptor
from ._graph_matrix import DistanceMatrix

__all__ = ("WienerIndex",)


class WienerIndex(Descriptor):
    r"""Wiener index.

    :type polarity: bool
    :param polarity: use polarity Wiener index
    """

    __slots__ = ("_polarity",)
    since = "1.0.0"
    explicit_hydrogens = False

    def description(self):
        return "Wiener {}index".format("polarity " if self._polarity else "")

    @classmethod
    def preset(cls, version):
        yield cls(False)
        yield cls(True)

    def __str__(self):
        return "WPol" if self._polarity else "WPath"

    def parameters(self):
        return (self._polarity,)

    def __init__(self, polarity=False):
        self._polarity = polarity

    def dependencies(self):
        return {"D": DistanceMatrix(self.explicit_hydrogens)}

    def calculate(self, D):
        if self._polarity:
            return int(0.5 * (D == 3).sum())
        else:
            return int(0.5 * D.sum())

    rtype = int
