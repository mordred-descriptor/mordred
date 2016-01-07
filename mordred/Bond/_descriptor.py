from ._bond_types import BondCount
from ._hbond import HBondAcceptor, HBondDonor
from ._rotatable import RotatableBondsCount, RotatableBondsRatio

_descriptors = [
    BondCount,
    HBondAcceptor, HBondDonor,
    RotatableBondsCount, RotatableBondsRatio
]

__all__ = [d.__name__ for d in _descriptors]
