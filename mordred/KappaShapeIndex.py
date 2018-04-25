from __future__ import division

from .Chi import ChiCache
from ._base import Descriptor

__all__ = ("KappaShapeIndex1", "KappaShapeIndex2", "KappaShapeIndex3")


class KappaShapeIndexBase(Descriptor):
    explicit_hydrogens = False

    __slots__ = ("_order",)

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return "Kier{}".format(self._order)

    def parameters(self):
        return ()

    def __init__(self):
        self._order = int(self.__class__.__name__[-1])

    def dependencies(self):
        return {"Chi": ChiCache(self._order)}

    def _common(self, Chi):
        P = len(Chi.path)

        A = self.mol.GetNumAtoms()
        Pmin = A - self._order

        return P, A, Pmin

    rtype = float


class KappaShapeIndex1(KappaShapeIndexBase):
    r"""Kappa shape index 1 descriptor.

    :returns: NaN when :math:`N_{\rm Chi-path} = 0`
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "kappa shape index 1"

    def calculate(self, Chi):
        P, A, Pmin = self._common(Chi)
        Pmax = 0.5 * A * (A - 1)

        with self.rethrow_zerodiv():
            return 2 * Pmax * Pmin / (P * P)


class KappaShapeIndex2(KappaShapeIndexBase):
    r"""Kappa shape index 2 descriptor.

    :returns: NaN when :math:`N_{\rm Chi-path} = 0`
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "kappa shape index 2"

    def calculate(self, Chi):
        P, A, Pmin = self._common(Chi)
        Pmax = 0.5 * (A - 1) * (A - 2)

        with self.rethrow_zerodiv():
            return 2 * Pmax * Pmin / (P * P)


class KappaShapeIndex3(KappaShapeIndexBase):
    r"""Kappa shape index 3 descriptor.

    :returns: NaN when :math:`N_{\rm Chi-path} = 0`
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "kappa shape index 3"

    def calculate(self, Chi):
        P, A, Pmin = self._common(Chi)

        if A % 2 == 0:
            Pmax = 0.25 * (A - 2) ** 2
        else:
            Pmax = 0.25 * (A - 1) * (A - 3)

        with self.rethrow_zerodiv():
            return 4 * Pmax * Pmin / (P * P)
