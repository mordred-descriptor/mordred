import numpy as np

from ._base import Descriptor
from ._common import AdjacencyMatrix


class WalkCount(Descriptor):
    r"""walk count descriptor.

    :type order: int
    :param order: walk length

    :type total: bool
    :param total: sum of walk count over 1 to order

    :type self_returning: bool
    :param self_returning: use self returning walk only

    :rtype: int
    """

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        for start, sr in [(1, False), (2, True)]:
            for l in range(start, 11):
                yield cls(l, False, sr)

            yield cls(10, True, sr)

    def __str__(self):
        T = '{}SRW{:02d}' if self._self_returning else '{}MWC{:02d}'
        return T.format('T' if self._total else '', self._order)

    __slots__ = ('_order', '_total', '_self_returning',)

    def __init__(self, order=1, total=False, self_returning=False):
        self._order = order
        self._total = total
        self._self_returning = self_returning

    def dependencies(self):
        if self._total:
            W = ('W', self.__class__(
                self._order,
                False,
                self._self_returning,
            ))

            if self._order > 1:
                T = ('T', self.__class__(
                    self._order - 1,
                    True,
                    self._self_returning,
                ))

                return dict([W, T])

            return dict([W])

        return dict(
            An=AdjacencyMatrix(
                self.explicit_hydrogens,
                order=self._order,
            )
        )

    def calculate(self, mol, An=None, T=None, W=None):
        if self._total:
            if self._order == 1:
                return mol.GetNumAtoms() + W

            return T + W

        if self._self_returning:
            return np.log(An.trace() + 1)

        else:
            if self._order == 1:
                return An.sum() / 2

            return np.log(An.sum() + 1)
