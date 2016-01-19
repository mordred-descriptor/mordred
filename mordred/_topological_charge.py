from itertools import chain

import numpy as np

from six import integer_types

from ._base import Descriptor
from ._common import AdjacencyMatrix, DistanceMatrix


class ChargeTermMatrix(Descriptor):
    explicit_hydrogens = False
    require_connected = False

    def dependencies(self):
        return dict(
            A=AdjacencyMatrix(
                self.explicit_hydrogens,
                False,
            ),
            D=DistanceMatrix(
                self.explicit_hydrogens,
                False,
                False,
            ),
        )

    def calculate(self, mol, A, D):
        D2 = D.copy()
        D2[D2 != 0] **= -2
        np.fill_diagonal(D2, 0)

        M = A.dot(D2)
        return M - M.T


class TopologicalCharge(Descriptor):
    r"""topological charge descriptor.

    :type type: str
    :param type:
        * 'sum': sum of order-distance atom pairs coefficient
        * 'mean': mean of order-distance atom pairs coefficient
        * 'global': sum of mean-topoCharge over 0 to order

    :type order: int
    :param order: int

    :rtype: float

    References
        * :cite:`10.1021/ci00019a008`
    """

    explicit_hydrogens = False
    require_connected = False

    @classmethod
    def preset(cls):
        return chain(
            (cls(t, o) for t in ('raw', 'mean') for o in range(1, 11)),
            [cls('global', 10)]
        )

    def __str__(self):
        if self.type == 'global':
            return 'JGT{}'.format(self.order)
        elif self.type == 'mean':
            return 'JGI{}'.format(self.order)
        else:
            return 'GGI{}'.format(self.order)

    descriptor_keys = 'type', 'order'

    def __init__(self, type='global', order=10):
        assert type in ['global', 'mean', 'raw']
        assert type == 'global' or isinstance(order, integer_types)

        self.type = type
        self.order = order

    def dependencies(self):
        return dict(
            CT=ChargeTermMatrix(),
            D=DistanceMatrix(
                self.explicit_hydrogens,
                False,
                False,
            )
        )

    def calculate(self, mol, CT, D):
        D = D * np.tri(*D.shape)
        D[D == 0] = np.inf

        f = D <= self.order if self.type == 'global' else D == self.order

        CT = CT[f]

        if self.type == 'raw':
            return np.abs(CT).sum()

        # create frequency vector
        Df = D[f]
        C = Df.copy()
        for i in np.unique(Df):
            C[Df == i] = len(Df[Df == i])

        return np.abs(CT / C).sum()
