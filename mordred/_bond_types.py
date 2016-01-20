from enum import IntEnum
from itertools import chain

from rdkit import Chem

from ._base import Descriptor, parse_enum


class BondType(IntEnum):
    any = 1
    heavy = 2

    single = 3
    double = 4
    triple = 5

    aromatic = 6
    multiple = 7


bond_types = (
    (BondType.any, ('', lambda _: True)),
    (BondType.heavy, ('O', lambda _: True)),

    (BondType.single, ('S', lambda b: b.GetBondType() == Chem.BondType.SINGLE)),
    (BondType.double, ('D', lambda b: b.GetBondType() == Chem.BondType.DOUBLE)),
    (BondType.triple, ('T', lambda b: b.GetBondType() == Chem.BondType.TRIPLE)),

    (BondType.aromatic, ('A', lambda b: b.GetIsAromatic() or
                         b.GetBondType() == Chem.BondType.AROMATIC)),
    (BondType.multiple, ('M', lambda b: b.GetIsAromatic() or
                         b.GetBondType() != Chem.BondType.SINGLE)),
)

bond_type_dict = dict(bond_types)


class BondCount(Descriptor):
    r"""bond count descriptor.

    :type type: str
    :param type: one of bond_types

    :type kekulize: bool
    :param kekulize: use kekulized structure

    :rtype: int
    """

    bond_types = tuple(b.name for b in BondType)

    require_connected = False

    @classmethod
    def preset(cls):
        return chain(
            map(lambda t: cls(t, False), BondType),
            map(lambda t: cls(t, True), [
                BondType.single, BondType.double,
            ])
        )

    def __str__(self):
        K = 'K' if self.kekulize else ''
        return 'nBonds{}{}'.format(K, self.bond_name)

    @property
    def explicit_hydrogens(self):
        return self.type in (BondType.any, BondType.single)

    descriptor_keys = 'type', 'kekulize'

    def __init__(self, type='any', kekulize=False):
        self.type = parse_enum(BondType, type)
        self.bond_name, self.check_bond = bond_type_dict[self.type]
        self.kekulize = kekulize

    def calculate(self, mol):
        return sum(1 for b in mol.GetBonds() if self.check_bond(b))
