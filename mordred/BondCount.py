from enum import IntEnum
from itertools import chain

from rdkit import Chem

from ._base import Descriptor
from ._util import parse_enum

__all__ = ("BondCount",)


class BondType(IntEnum):
    __slots__ = ()

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
    (BondType.any, ("", "all bonds", lambda _: True)),
    (BondType.heavy, ("O", "bonds connecting to heavy atom", lambda _: True)),
    (
        BondType.single,
        ("S", "single bonds", lambda b: b.GetBondType() == Chem.BondType.SINGLE),
    ),
    (
        BondType.double,
        ("D", "double bonds", lambda b: b.GetBondType() == Chem.BondType.DOUBLE),
    ),
    (
        BondType.triple,
        ("T", "triple bonds", lambda b: b.GetBondType() == Chem.BondType.TRIPLE),
    ),
    (
        BondType.aromatic,
        (
            "A",
            "aromatic bonds",
            lambda b: b.GetIsAromatic() or b.GetBondType() == Chem.BondType.AROMATIC,
        ),
    ),
    (
        BondType.multiple,
        (
            "M",
            "multiple bonds",
            lambda b: b.GetIsAromatic() or b.GetBondType() != Chem.BondType.SINGLE,
        ),
    ),
)

bond_type_dict = dict(bond_types)


class BondCount(Descriptor):
    r"""bond count descriptor.

    :type type: str
    :param type: one of bond_types

    :type kekulize: bool
    :param kekulize: use kekulized structure
    """

    since = "1.0.0"
    __slots__ = ("_type", "_bond_name", "_bond_desc", "_check_bond", "kekulize")

    def description(self):
        return "number of {} in {}kekulized structure".format(
            self._bond_desc, "" if self.kekulize else "non-"
        )

    bond_types = tuple(b.name for b in BondType)

    @classmethod
    def preset(cls, version):
        return chain(
            map(lambda t: cls(t, False), BondType),
            map(lambda t: cls(t, True), [BondType.single, BondType.double]),
        )

    def __str__(self):
        K = "K" if self.kekulize else ""
        return "nBonds{}{}".format(K, self._bond_name)

    @property
    def explicit_hydrogens(self):
        return self._type in (BondType.any, BondType.single)

    def parameters(self):
        return self._type, self.kekulize

    def __init__(self, type="any", kekulize=False):
        self._type = parse_enum(BondType, type)
        self._bond_name, self._bond_desc, self._check_bond = bond_type_dict[self._type]
        self.kekulize = kekulize

    def calculate(self):
        return sum(1 for b in self.mol.GetBonds() if self._check_bond(b))

    rtype = int

    _extra_docs = ("bond_types",)
