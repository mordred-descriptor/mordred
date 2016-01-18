from numpy import nan

from ._base import Descriptor
from ._chi import ChiCache


class KappaShapeIndex(Descriptor):
    r"""kappa shape index descriptor.

    :type order: int
    :param order: order of kier, [1,3]

    :rtype: float
    """

    explicit_hydrogens = False
    require_connected = False

    @classmethod
    def preset(cls):
        return map(cls, range(1, 4))

    def __str__(self):
        return 'Kier{}'.format(self.order)

    descriptor_keys = 'order',

    def __init__(self, order=1):
        assert order in [1, 2, 3]

        self.order = order

    def dependencies(self):
        return dict(
            Chi=ChiCache(self.order)
        )

    def calculate(self, mol, Chi):
        P = len(Chi.path)
        if P == 0:
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
