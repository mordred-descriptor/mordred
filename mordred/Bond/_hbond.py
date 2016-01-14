from .._base import Descriptor
from rdkit.Chem import rdMolDescriptors


class HBondAcceptor(Descriptor):
    '''
    hydrogen bond acceptor descriptor

    Returns:
        int: hydrogen bond acceptor count
    '''

    explicit_hydrogens = False
    descriptor_name = 'nHBAcc'

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBA(mol)


class HBondDonor(Descriptor):
    '''
    hydrogen bond donor descriptor

    Returns:
        int: hydrogen bond donor count
    '''

    explicit_hydrogens = False
    descriptor_name = 'nHBDon'

    def calculate(self, mol):
        return rdMolDescriptors.CalcNumHBD(mol)
