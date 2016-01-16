from ._base import Descriptor
from rdkit.Chem import rdMolDescriptors


class HBondAcceptor(Descriptor):
    r'''
    hydrogen bond acceptor descriptor

    Returns:
        int: hydrogen bond acceptor count
    '''

    explicit_hydrogens = False

    def __str__(self):
        return 'nHBAcc'

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBA(mol)


class HBondDonor(Descriptor):
    r'''
    hydrogen bond donor descriptor

    Returns:
        int: hydrogen bond donor count
    '''

    explicit_hydrogens = False

    def __str__(self):
        return 'nHBDon'

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBD(mol)
