from ._base import Descriptor
from rdkit.Chem import BondType
from itertools import chain


bond_type_dict = {
    'heavy': ('O', lambda _: True),
    'any': ('', lambda _: True),

    'single': ('S', lambda b: b.GetBondType() == BondType.SINGLE),
    'double': ('D', lambda b: b.GetBondType() == BondType.DOUBLE),
    'triple': ('T', lambda b: b.GetBondType() == BondType.TRIPLE),
    'aromatic': ('A', lambda b: b.GetIsAromatic() or b.GetBondType() == BondType.AROMATIC),

    'multiple': ('M', lambda b: b.GetIsAromatic() or b.GetBondType() != BondType.SINGLE),
}


class BondCount(Descriptor):
    '''
    bond count descriptor

    Parameters:
        bond_type(str):
            * 'any' - any
            * 'heavy' - any, not include bond which connect to hydrogen

            * 'single' - single bonds
            * 'double' - double bonds
            * 'triple' - triple bonds
            * 'aromatic' - aromatic bonds

            * 'multiple' - multiple, include aromatic

        kekulize(bool): use kekulized structure

    Returns:
        int: bond count
    '''

    @classmethod
    def preset(cls):
        return chain(
            map(lambda t: cls(t, False), [
                'any', 'heavy',
                'single', 'double', 'triple', 'aromatic',
                'multiple',
            ]),
            map(lambda t: cls(t, True), [
                'single', 'double',
            ])
        )

    @property
    def descriptor_name(self):
        K = 'K' if self.kekulize else ''
        return 'nBonds{}{}'.format(K, self.bond_name)

    @property
    def explicit_hydrogens(self):
        return self.bond_name in ['', 'S']

    @property
    def descriptor_key(self):
        return self.make_key(self.bond_name, self.kekulize)

    def __init__(self, bond_type='any', kekulize=False):
        self.bond_name, self.check_bond = bond_type_dict[bond_type.lower()]
        self.kekulize = kekulize

    def calculate(self, mol):
        return sum(1 for b in mol.GetBonds() if self.check_bond(b))
