from .._base import Descriptor
from ..Chi._descriptor import ChiCache
from numpy import nan


class KappaShapeIndex(Descriptor):
    explicit_hydrogens = False

    descriptor_defaults = [(1,), (2,), (3,)]

    def __init__(self, order=1):
        assert order in [1,2,3]

        self.order = order

    @property
    def dependencies(self):
        return dict(
            Chi=ChiCache.make_key(self.order)
        )

    @property
    def descriptor_name(self):
        return 'Kier{}'.format(self.order)

    @property
    def descriptor_key(self):
        return self.make_key(self.order)

    def calculate(self, mol, Chi):
        P = len(Chi.path)
        if P <= 0:
            return nan

        A = mol.GetNumAtoms()
        Pmin = A - self.order
        if self.order == 1:
            Pmax = float(A * (A - 1)) / 2.0
            return 2 * Pmax * Pmin / (P * P)

        elif self.order == 2:
            Pmax = float((A - 1) * (A - 2)) / 2.0
            return 2 * Pmax * Pmin / (P * P)

        elif A % 2 == 0:
            Pmax = float((A - 2) ** 2) / 4.0
        else:
            Pmax = float((A - 1) * (A - 3)) / 4.0

        return 4 * Pmax * Pmin / (P * P)

_descriptors = [KappaShapeIndex]
__all__ = [d.__name__ for d in _descriptors]
