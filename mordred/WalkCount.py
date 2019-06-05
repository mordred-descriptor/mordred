import numpy as np

from ._base import Descriptor
from ._graph_matrix import AdjacencyMatrix

__all__ = ("WalkCount",)


class WalkCount(Descriptor):
    r"""walk count descriptor.

    :type order: int
    :param order: walk length

    :type total: bool
    :param total: sum of walk count over 1 to order

    :type self_returning: bool
    :param self_returning: use self returning walk only
    """

    since = "1.0.0"
    __slots__ = ("_order", "_total", "_self_returning")

    explicit_hydrogens = False

    def description(self):
        return "{}walk count (leg-{}{})".format(
            "total " if self._total else "",
            self._order,
            ", only self returning walk" if self._self_returning else "",
        )

    @classmethod
    def preset(cls, version):
        for start, sr in [(1, False), (2, True)]:
            for l in range(start, 11):
                yield cls(l, False, sr)

            yield cls(10, True, sr)

    def __str__(self):
        T = "{}SRW{:02d}" if self._self_returning else "{}MWC{:02d}"
        return T.format("T" if self._total else "", self._order)

    def parameters(self):
        return self._order, self._total, self._self_returning

    def __init__(self, order=1, total=False, self_returning=False):
        self._order = order
        self._total = total
        self._self_returning = self_returning

    def dependencies(self):
        if self._total:
            d = {}
            d["W"] = self.__class__(self._order, False, self._self_returning)

            if self._order > 1:
                d["T"] = self.__class__(self._order - 1, True, self._self_returning)

            return d

        return {"An": AdjacencyMatrix(self.explicit_hydrogens, order=self._order)}

    def calculate(self, An=None, T=None, W=None):
        if self._total:
            if self._order == 1:
                return self.mol.GetNumAtoms() + W

            return T + W

        if self._self_returning:
            return np.log(An.trace() + 1)

        else:
            if self._order == 1:
                return 0.5 * An.sum()

            return np.log(An.sum() + 1)

    rtype = float
