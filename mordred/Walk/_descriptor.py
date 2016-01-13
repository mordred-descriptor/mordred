from .._base import Descriptor
from .._common import AdjacencyMatrix
import numpy as np


class WalkCount(Descriptor):
    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        for start, sr in [(1, False), (2, True)]:
            for l in range(start, 11):
                yield l, False, sr

            yield 10, True, sr

    def __init__(self, order=2, total=False, self_returning=False):
        self.order = order
        self.total = total
        self.self_returning = self_returning

    @property
    def dependencies(self):
        if self.total:
            W = ('W', self.make_key(
                self.order,
                False,
                self.self_returning,
            ))

            if self.order > 1:
                T = ('T', self.make_key(
                    self.order - 1,
                    True,
                    self.self_returning
                ))

                return dict([W, T])

            return dict([W])

        return dict(
            An=AdjacencyMatrix.make_key(
                self.explicit_hydrogens,
                False,
                self.order,
            )
        )

    @property
    def descriptor_name(self):
        T = '{}SRW{:02d}' if self.self_returning else '{}MWC{:02d}'
        return T.format('T' if self.total else '', self.order)

    @property
    def descriptor_key(self):
        return self.make_key(
            self.order,
            self.total,
            self.self_returning,
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

_descriptors = [WalkCount]
__all__ = [d.__name__ for d in _descriptors]
