from ._base import Descriptor
from ._common import Valence
import numpy as np


class ZagrebIndex(Descriptor):
    r'''
    Zagreb index descriptor

    .. math::

        {}^\lambda M_1 = \sum_{atoms} d_i^\lambda

        {}^\lambda M_2 = \sum_{edges} \left(d_i \cdot d_j \right)^\lambda

    where
    :math:`d_i` is degree of i-th atom

    :type version: int
    :param version: Zagreb index version. 1 or 2.

    :type variable: int
    :param variable: lambda value.

    :rtype: int
    '''

    explicit_hydrogens = False
    require_connected = False

    @classmethod
    def preset(cls):
        return (cls(v, x) for x in [1, -1] for v in [1, 2])

    def __str__(self):
        if self.variable in [1, -1]:
            m = '' if self.variable == 1 else 'm'
            return '{}Zagreb{}'.format(m, self.version)

        return 'Zagreb{}_{}'.format(self.version, self.variable)

    descriptor_keys = 'version', 'variable'

    def __init__(self, version=1, variable=1):
        assert version in [1, 2]
        self.version = version
        self.variable = variable

    def dependencies(self):
        return dict(
            V=Valence(
                self.explicit_hydrogens,
                False,
            ),
        )

    def calculate(self, mol, V):
        if not isinstance(self.variable, int) or self.variable < 0:
            V = V.astype('float')

        if self.version == 1:
            if np.any(V == 0):
                return np.nan

            return (V ** (self.variable * 2)).sum()
        else:
            return sum(
                (V[b.GetBeginAtomIdx()] * V[b.GetEndAtomIdx()]) ** self.variable
                for b in mol.GetBonds()
            )
