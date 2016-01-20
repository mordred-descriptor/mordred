from ._base import Descriptor
from ._common import DistanceMatrix


class WienerIndex(Descriptor):
    r"""Wiener index.

    :type polarity: bool
    :param polarity: use polarity Wiener index

    :rtype: int
    """

    explicit_hydrogens = False
    require_connected = False

    @classmethod
    def preset(cls):
        yield cls(False)
        yield cls(True)

    def __str__(self):
        return 'WPol' if self.polarity else 'WPath'

    __slots__ = ('polarity',)

    def __init__(self, polarity=False):
        self.polarity = polarity

    def dependencies(self):
        return dict(
            D=DistanceMatrix(
                self.explicit_hydrogens,
                False,
                False,
            )
        )

    def calculate(self, mol, D):
        if self.polarity:
            return int(0.5 * (D == 3).sum())
        else:
            return int(0.5 * D.sum())
