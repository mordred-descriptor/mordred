from ._base import Descriptor
from ._common import AdjacencyMatrix
import numpy as np


class WalkCount(Descriptor):
    r'''
    walk count descriptor

    Parameters:
        order(int): walk length
        total(bool): sum of walk count over 1 to order
        self_returning(bool): use self returning walk only

    Returns:
        int: walk count
    '''

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        for start, sr in [(1, False), (2, True)]:
            for l in range(start, 11):
                yield cls(l, False, sr)

            yield cls(10, True, sr)

    def __str__(self):
        T = '{}SRW{:02d}' if self.self_returning else '{}MWC{:02d}'
        return T.format('T' if self.total else '', self.order)

    descriptor_keys = 'order', 'total', 'self_returning'

    def __init__(self, order=1, total=False, self_returning=False):
        self.order = order
        self.total = total
        self.self_returning = self_returning

    @property
    def dependencies(self):
        if self.total:
            W = ('W', self.__class__(
                self.order,
                False,
                self.self_returning,
            ))

            if self.order > 1:
                T = ('T', self.__class__(
                    self.order - 1,
                    True,
                    self.self_returning,
                ))

                return dict([W, T])

            return dict([W])

        return dict(
            An=AdjacencyMatrix(
                self.explicit_hydrogens,
                False,
                self.order,
            )
        )

    def calculate(self, mol, An=None, T=None, W=None):
        if self.total:
            if self.order == 1:
                return mol.GetNumAtoms() + W

            return T + W

        if self.self_returning:
            return np.log(An.trace() + 1)

        else:
            if self.order == 1:
                return An.sum() / 2

            return np.log(An.sum() + 1)
