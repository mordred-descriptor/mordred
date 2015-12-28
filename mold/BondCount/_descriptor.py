from .._base import Descriptor


class BondCount(Descriptor):
    descriptor_defaults = [
        ('',), ('2',),
        ('S',), ('S2',), ('S3',),
        ('D',), ('D2',),
        ('T',), ('Q',), ('M',),
    ]

    @property
    def descriptor_name(self):
        return 'nBonds{}'.format(self.bond_type)

    @property
    def explicitHydrogens(self):
        return self.bond_type in set(['2', 'S', 'S2'])

    @property
    def kekulize(self):
        return self.bond_type != 'M'

    @property
    def descriptor_key(self):
        return self.make_key(self.bond_type)

    def __init__(self, bond_type):
        self.bond_type = bond_type

    def nBonds(self, n, mol):
        return sum((1 for b in mol.GetBonds()
                    if b.GetBondTypeAsDouble() == n))

    def nBondsNA(self, n, mol):
        return sum((1 for b in mol.GetBonds()
                    if b.GetBondTypeAsDouble() == n
                    and not b.GetIsAromatic()))

    def nBondsM(self, mol):
        return sum((1 for b in mol.GetBonds()
                    if b.GetBondTypeAsDouble() > 1.0))

    def calculate(self, mol):
        bt = self.bond_type
        if bt in ' 2':
            return mol.GetNumBonds()
        elif bt == 'S':
            return self.nBonds(1, mol)
        elif bt in ['S2', 'S3']:
            return self.nBondsNA(1, mol)
        elif bt == 'D':
            return self.nBonds(2, mol)
        elif bt == 'D2':
            return self.nBondsNA(2, mol)
        elif bt == 'T':
            return self.nBonds(3, mol)
        elif bt == 'Q':
            return self.nBonds(4, mol)
        elif bt == 'M':
            return self.nBondsM(mol)

_descriptors = [BondCount]
__all__ = [d.__name__ for d in _descriptors]
