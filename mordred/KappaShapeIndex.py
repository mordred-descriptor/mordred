from numpy import nan

from ._base import Descriptor
from .Chi import ChiCache


__all__ = ('KappaShapeIndex1', 'KappaShapeIndex2', 'KappaShapeIndex3',)


class KappaShapeIndexBase(Descriptor):
    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return 'Kier{}'.format(self._order)

    def as_key(self):
        return self.__class__, ()

    def __init__(self):
        self._order = int(self.__class__.__name__[-1])

    def dependencies(self):
        return {'Chi': ChiCache(self._order)}

    def _common(self, mol, Chi):
        P = len(Chi.path)
        if P == 0:
            P = nan

        A = mol.GetNumAtoms()
        Pmin = A - self._order

        return P, A, Pmin

    rtype = float


class KappaShapeIndex1(KappaShapeIndexBase):
    r"""Kappa shape index 1 descriptor.

    :returns: NaN when :math:`N_{\rm Chi-path} = 0`
    """

    __slots__ = ()

    def calculate(self, mol, Chi):
        P, A, Pmin = self._common(mol, Chi)
        Pmax = float(A * (A - 1)) / 2.0

        return 2 * Pmax * Pmin / (P * P)


class KappaShapeIndex2(KappaShapeIndexBase):
    r"""Kappa shape index 2 descriptor.

    :returns: NaN when :math:`N_{\rm Chi-path} = 0`
    """

    __slots__ = ()

    def calculate(self, mol, Chi):
        P, A, Pmin = self._common(mol, Chi)
        Pmax = float((A - 1) * (A - 2)) / 2.0

        return 2 * Pmax * Pmin / (P * P)


class KappaShapeIndex3(KappaShapeIndexBase):
    r"""Kappa shape index 3 descriptor.

    :returns: NaN when :math:`N_{\rm Chi-path} = 0`
    """

    __slots__ = ()

    def calculate(self, mol, Chi):
        P, A, Pmin = self._common(mol, Chi)

        if A % 2 == 0:
            Pmax = float((A - 2) ** 2) / 4.0
        else:
            Pmax = float((A - 1) * (A - 3)) / 4.0

        return 4 * Pmax * Pmin / (P * P)
