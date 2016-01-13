from .._base import Descriptor


class AromaticAtomsCount(Descriptor):
    descriptor_name = 'nAromAtom'

    def calculate(self, mol):
        return sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())


class AromaticBondsCount(Descriptor):
    descriptor_name = 'nAromBond'

    def calculate(self, mol):
        return sum(1 for b in mol.GetBonds() if b.GetIsAromatic())

_descriptors = [AromaticAtomsCount, AromaticBondsCount]
__all__ = [d.__name__ for d in _descriptors]
