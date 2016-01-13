from .._base import Descriptor
from rdkit.Chem import rdMolDescriptors


class HBondAcceptor(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'nHBAcc'

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBA(mol)


class HBondDonor(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'nHBDon'

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBD(mol)
