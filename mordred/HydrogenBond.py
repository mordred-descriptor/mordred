from rdkit.Chem import rdMolDescriptors

from ._base import Descriptor

__all__ = ("HBondAcceptor", "HBondDonor")


class HBondBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False

    @classmethod
    def preset(cls, version):
        yield cls()

    rtype = int


class HBondAcceptor(HBondBase):
    r"""hydrogen bond acceptor descriptor(rdkit wrapper)."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "number of hydrogen bond acceptor"

    def __str__(self):
        return "nHBAcc"

    def parameters(self):
        return ()

    def calculate(self):
        return rdMolDescriptors.CalcNumHBA(self.mol)


class HBondDonor(HBondBase):
    r"""hydrogen bond donor descriptor(rdkit wrapper)."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "number of hydrogen bond donor"

    def __str__(self):
        return "nHBDon"

    def parameters(self):
        return ()

    def calculate(self):
        return rdMolDescriptors.CalcNumHBD(self.mol)
