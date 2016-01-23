from collections import defaultdict

import numpy as np

from ._base import Descriptor


class CarbonTypesBase(Descriptor):
    explicit_hydrogens = False
    kekulize = True
    require_connected = False


class CarbonTypesCache(CarbonTypesBase):
    __slots__ = ()

    def calculate(self, mol):
        r = defaultdict(lambda: defaultdict(int))
        for a in mol.GetAtoms():
            if a.GetAtomicNum() != 6:
                continue

            double = 0
            triple = 0
            carbon = 0

            for b in a.GetBonds():
                other = b.GetBeginAtom()
                if a.GetIdx() == other.GetIdx():
                    other = b.GetEndAtom()

                if other.GetAtomicNum() == 6:
                    carbon += 1

                bt = b.GetBondTypeAsDouble()
                if bt == 2.0:
                    double += 1
                elif bt == 3.0:
                    triple += 1

            if (double == 2 and triple == 0) or (triple == 1 and double == 0):
                SP = 1
            elif double == 1 and triple == 0:
                SP = 2
            else:
                SP = 3

            r[SP][carbon] += 1

        return r


class CarbonTypes(CarbonTypesBase):
    r"""carbon types descriptor.

    :type nCarbon: int
    :param nCarbon: count `n`-carbon bonded carbon

    :type SP: int
    :param SP: count :math:`{\rm SP}n` carbon

    :rtype: int
    """

    @classmethod
    def preset(cls):
        return map(lambda args: cls(*args), [
            (1, 1), (2, 1),
            (1, 2), (2, 2), (3, 2),
            (1, 3), (2, 3), (3, 3), (4, 3),
        ])

    def __str__(self):
        return 'C{}SP{}'.format(self._nCarbon, self._SP)

    __slots__ = ('_nCarbon', '_SP',)

    def __init__(self, nCarbon=1, SP=3):
        assert SP in [1, 2, 3]

        self._nCarbon = nCarbon
        self._SP = SP

    def dependencies(self):
        return dict(CT=CarbonTypesCache())

    def calculate(self, mol, CT):
        return CT[self._SP][self._nCarbon]


class HybridizationRatio(CarbonTypesBase):
    r"""hybridization ratio descriptor.

    .. math::

        {\rm HybRatio} = \frac{N_{\rm SP3}}{N_{\rm SP2} + N_{\rm SP3}}

    :rtype: float
    :returns: NaN when :math:`N_{\rm SP2} + N_{\rm SP3} = 0`.
    """

    __slots__ = ()

    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return 'HybRatio'

    def dependencies(self):
        return dict(CT=CarbonTypesCache())

    def calculate(self, mol, CT):
        Nsp3 = float(sum(CT[3].values()))
        Nsp2 = float(sum(CT[2].values()))

        if Nsp3 == Nsp2 == 0:
            return np.nan

        return Nsp3 / (Nsp2 + Nsp3)
