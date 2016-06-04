from enum import IntEnum
from itertools import chain

from rdkit import Chem

from ._base import Descriptor
from ._util import parse_enum


__all__ = ('BondCount',)


class BondType(IntEnum):
    any = 1
    heavy = 2

    single = 3
    double = 4
    triple = 5

    aromatic = 6
    multiple = 7

    @property
    def as_argument(self):
        return self.name


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
    """

    bond_types = tuple(b.name for b in BondType)

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
        return 'nBonds{}{}'.format(K, self._bond_name)

    @property
    def explicit_hydrogens(self):
        return self._type in (BondType.any, BondType.single)

    __slots__ = ('_type', '_bond_name', '_check_bond', 'kekulize',)

    def as_key(self):
        return self.__class__, (self._type, self.kekulize)

    def __init__(self, type='any', kekulize=False):
        self._type = parse_enum(BondType, type)
        self._bond_name, self._check_bond = bond_type_dict[self._type]
        self.kekulize = kekulize

    def calculate(self, mol):
        return sum(1 for b in mol.GetBonds() if self._check_bond(b))

    rtype = int
