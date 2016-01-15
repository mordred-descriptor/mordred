from ._base import Descriptor


class AromaticAtomsCount(Descriptor):
    '''
    aromatic atoms count descriptor

    Returns:
        int: count of aromatic atoms
    '''

    descriptor_name = 'nAromAtom'

    def calculate(self, mol):
        return sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())


class AromaticBondsCount(Descriptor):
    '''
    aromatic bonds count descriptor

    Returns:
        int: count of aromatic bonds
    '''

    descriptor_name = 'nAromBond'

    def calculate(self, mol):
        return sum(1 for b in mol.GetBonds() if b.GetIsAromatic())

_descriptors = [AromaticAtomsCount, AromaticBondsCount]
__all__ = [d.__name__ for d in _descriptors]
