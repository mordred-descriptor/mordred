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
    r'''
    bond count descriptor

    :type type: str
    :param type:
        * 'any' - any
        * 'heavy' - any, not include bond which connect to hydrogen

        * 'single' - single bonds
        * 'double' - double bonds
        * 'triple' - triple bonds
        * 'aromatic' - aromatic bonds

        * 'multiple' - multiple, include aromatic

    :type kekulize: bool
    :param kekulize: use kekulized structure

    :rtype: int
    '''

    require_connected = False

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

    def __str__(self):
        K = 'K' if self.kekulize else ''
        return 'nBonds{}{}'.format(K, self.bond_name)

    @property
    def explicit_hydrogens(self):
        return self.bond_name in ['', 'S']

    descriptor_keys = 'type', 'kekulize'

    def __init__(self, type='any', kekulize=False):
        self.type = type
        self.bond_name, self.check_bond = bond_type_dict[type.lower()]
        self.kekulize = kekulize

    def calculate(self, mol):
        return sum(1 for b in mol.GetBonds() if self.check_bond(b))
