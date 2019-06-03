from ._base import Descriptor
from .SLogP import SMR, SLogP
from .Weight import Weight
from .HydrogenBond import HBondDonor, HBondAcceptor

__all__ = ("Lipinski", "GhoseFilter")


class LipinskiLike(Descriptor):
    __slots__ = ()

    @classmethod
    def preset(cls, version):
        yield cls()

    def parameters(self):
        return ()

    def __str__(self):
        return self.__class__.__name__

    rtype = bool


class Lipinski(LipinskiLike):
    r"""Lipinski rule of 5 descriptor.

    LogP: SLogP
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "Lipinski rule of five"

    def dependencies(self):
        return {
            "HBAcc": HBondAcceptor(),
            "HBDon": HBondDonor(),
            "LogP": SLogP(),
            "MW": Weight(),
        }

    def calculate(self, LogP, MW, HBDon, HBAcc):
        return HBDon <= 5 and HBAcc <= 10 and MW <= 500 and LogP <= 5


class GhoseFilter(LipinskiLike):
    r"""Ghose filter descriptor.

    LogP, MR: SLogP, SMR
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "Ghose filter"

    def dependencies(self):
        return {"LogP": SLogP(), "MR": SMR(), "MW": Weight()}

    def calculate(self, MW, LogP, MR):
        return (
            (160 <= MW <= 480)
            and (20 <= self.mol.GetNumAtoms() <= 70)
            and (-0.4 <= LogP <= 5.6)
            and (40 <= MR <= 130)
        )
