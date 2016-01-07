from .._base import Descriptor
from rdkit.Chem import rdMolDescriptors


class HBondAcceptor(Descriptor):
    explicit_hydrogens = False

    descriptor_defaults = [()]

    @property
    def descriptor_name(self):
        return 'nHBAcc'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBA(mol)


class HBondDonor(Descriptor):
    explicit_hydrogens = False

    descriptor_defaults = [()]

    @property
    def descriptor_name(self):
        return 'nHBDon'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBD(mol)
