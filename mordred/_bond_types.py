from ._base import Descriptor
from rdkit.Chem import BondType


class BondCount(Descriptor):
    '''
    bond count descriptor

    Parameters:
        bond_type(str):
            * '' - any
            * 'O' - any, not include bond which connect to hydrogen
            * 'T' - triple,
            * 'M' - multiple, include aromatic

    Returns:
        int: bond count
    '''

    @classmethod
    def preset(cls):
        return map(cls, [
            '', 'O', 'T', 'M',
        ])

    @property
    def descriptor_name(self):
        return 'nBonds{}'.format(self.bond_type)

    @property
    def explicit_hydrogens(self):
        return self.bond_type == ''

    @property
    def descriptor_key(self):
        return self.make_key(self.bond_type)

    def __init__(self, bond_type=''):
        assert bond_type in set(['', 'O', 'T', 'M'])
        self.bond_type = bond_type

    def calculate(self, mol):
        bt = self.bond_type
        if bt in ['', 'O']:
            return sum(1 for _ in mol.GetBonds())
        elif bt == 'M':
            return sum(1 for b in mol.GetBonds()
                       if b.GetBondTypeAsDouble() > 1 or b.GetIsAromatic())
        else:
            return sum(1 for b in mol.GetBonds() if b.GetBondType() == BondType.TRIPLE)
