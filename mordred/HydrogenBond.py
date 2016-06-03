from rdkit.Chem import rdMolDescriptors

from ._base import Descriptor


__all__ = (
    'HBondAcceptor', 'HBondDonor',
)


class HBondBase(Descriptor):
    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    rtype = int


class HBondAcceptor(HBondBase):
    r"""hydrogen bond acceptor descriptor(rdkit wrapper)."""

    __slots__ = ()

    def __str__(self):
        return 'nHBAcc'

    def as_key(self):
        return self.__class__, ()

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBA(mol)


class HBondDonor(HBondBase):
    r"""hydrogen bond donor descriptor(rdkit wrapper)."""

    __slots__ = ()

    def __str__(self):
        return 'nHBDon'

    def as_key(self):
        return self.__class__, ()

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBD(mol)
