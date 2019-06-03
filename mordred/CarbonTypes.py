from __future__ import division

from collections import defaultdict

from rdkit.Chem import HybridizationType
from rdkit.Chem.rdMolDescriptors import CalcFractionCSP3

from ._base import Descriptor

__all__ = ("CarbonTypes", "HybridizationRatio", "FractionCSP3")


class CarbonTypesBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False
    kekulize = True


class CarbonTypesCache(CarbonTypesBase):
    __slots__ = ()

    def parameters(self):
        return ()

    _hybridization = {
        HybridizationType.SP: 1,
        HybridizationType.SP2: 2,
        HybridizationType.SP3: 3,
        HybridizationType.SP3D: 3,
        HybridizationType.SP3D2: 3,
    }

    def calculate(self):
        r = defaultdict(lambda: defaultdict(int))
        for a in self.mol.GetAtoms():
            if a.GetAtomicNum() != 6:
                continue

            carbon = sum(other.GetAtomicNum() == 6 for other in a.GetNeighbors())

            SP = self._hybridization.get(a.GetHybridization())

            r[SP][carbon] += 1

        return r


class CarbonTypes(CarbonTypesBase):
    r"""carbon types descriptor.

    :type nCarbon: int
    :param nCarbon: count `n`-carbon bonded carbon

    :type SP: int
    :param SP: count :math:`{\rm SP}n` carbon
    """

    since = "1.0.0"
    __slots__ = ("_nCarbon", "_SP")

    def description(self):
        return "SP{} carbon bound to {} other carbon{}".format(
            self._SP if self._SP != 1 else "",
            self._nCarbon,
            "s" if self._nCarbon > 1 else "",
        )

    @classmethod
    def preset(cls, version):
        return map(
            lambda args: cls(*args),
            [(1, 1), (2, 1), (1, 2), (2, 2), (3, 2), (1, 3), (2, 3), (3, 3), (4, 3)],
        )

    def __str__(self):
        return "C{}SP{}".format(self._nCarbon, self._SP)

    def parameters(self):
        return self._nCarbon, self._SP

    def __init__(self, nCarbon=1, SP=3):
        assert SP in [1, 2, 3]

        self._nCarbon = nCarbon
        self._SP = SP

    def dependencies(self):
        return {"CT": CarbonTypesCache()}

    def calculate(self, CT):
        return CT[self._SP][self._nCarbon]

    rtype = int


class HybridizationRatio(CarbonTypesBase):
    r"""hybridization ratio descriptor.

    .. math::

        {\rm HybRatio} = \frac{N_{\rm SP3}}{N_{\rm SP2} + N_{\rm SP3}}

    :returns: NaN when :math:`N_{\rm SP2} + N_{\rm SP3} = 0`.
    """

    def description(self):
        return "hybridization ratio"

    since = "1.0.0"
    __slots__ = ()

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return "HybRatio"

    def parameters(self):
        return ()

    def dependencies(self):
        return {"CT": CarbonTypesCache()}

    def calculate(self, CT):
        Nsp3 = sum(CT[3].values())
        Nsp2 = sum(CT[2].values())

        if Nsp3 == Nsp2 == 0:
            self.fail(ValueError("there are no sp3 and sp2 carbons"))

        return Nsp3 / (Nsp2 + Nsp3)

    rtype = float


class FractionCSP3(Descriptor):
    r"""the fraction of C atoms that are SP3 hybridized."""

    __slots__ = ()
    since = "1.1.0"

    @classmethod
    def preset(cls, version):
        yield cls()

    def description(self):
        return "the fraction of C atoms that are SP3 hybridized"

    def __str__(self):
        return "FCSP3"

    def parameters(self):
        return ()

    def calculate(self):
        return CalcFractionCSP3(self.mol)

    rtype = float
