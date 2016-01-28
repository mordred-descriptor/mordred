import numpy as np

from ._base import Descriptor
from ._common import Valence


class ZagrebIndex(Descriptor):
    r"""Zagreb index descriptor.

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
    :returns: NaN when valence of any atoms are 0
    """

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        return (cls(v, x) for x in [1, -1] for v in [1, 2])

    def __str__(self):
        if self._variable in [1, -1]:
            m = '' if self._variable == 1 else 'm'
            return '{}Zagreb{}'.format(m, self._version)

        return 'Zagreb{}_{}'.format(self._version, self._variable)

    __slots__ = ('_version', '_variable',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._version, self._variable)

    def __init__(self, version=1, variable=1):
        assert version in [1, 2]
        self._version = version
        self._variable = variable

    def dependencies(self):
        return dict(
            V=Valence(self.explicit_hydrogens),
        )

    def calculate(self, mol, V):
        if not isinstance(self._variable, int) or self._variable < 0:
            V = V.astype('float')

        if self._version == 1:
            if np.any(V == 0):
                return np.nan

            return (V ** (self._variable * 2)).sum()
        else:
            return sum(
                (V[b.GetBeginAtomIdx()] * V[b.GetEndAtomIdx()]) ** self._variable
                for b in mol.GetBonds()
            )
