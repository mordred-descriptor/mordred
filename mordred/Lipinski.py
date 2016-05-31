from .HydrogenBond import HBondAcceptor, HBondDonor
from .SLogP import SLogP, SMR
from .Weight import Weight

from ._base import Descriptor


__all__ = ('Lipinski', 'GhoseFilter')


class LipinskiLike(Descriptor):

    @classmethod
    def preset(cls):
        yield cls()

    def __reduce_ex__(self, version):
        return self.__class__, ()

    def __str__(self):
        return self.__class__.__name__

    rtype = bool


class Lipinski(LipinskiLike):
    r"""Lipinski rule of 5 descriptor.

    LogP: SLogP
    """

    __slots__ = ()

    def dependencies(self):
        return {
            'LogP': SLogP(),
            'MW': Weight(),
            'HBDon': HBondDonor(),
            'HBAcc': HBondAcceptor(),
        }

    def calculate(self, mol, LogP, MW, HBDon, HBAcc):
        return\
            HBDon <= 5 and\
            HBAcc <= 10 and\
            MW <= 500 and\
            LogP <= 5


class GhoseFilter(LipinskiLike):
    r"""Ghose filter descriptor.

    LogP, MR: SLogP, SMR
    """

    __slots__ = ()

    def dependencies(self):
        return {
            'LogP': SLogP(),
            'MR': SMR(),
            'MW': Weight(),
        }

    def calculate(self, mol, MW, LogP, MR):
        return\
            (160 <= MW <= 480) and\
            (20 <= mol.GetNumAtoms() <= 70) and\
            (-0.4 <= LogP <= 5.6) and\
            (40 <= MR <= 130)
