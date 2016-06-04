from math import pi

import numpy as np

from .BondCount import BondCount
from .RingCount import RingCount

from ._base import Descriptor
from ._atomic_property import table


__all__ = (
    'VdwVolumeABC',
)


bondi_radii = {
    'H': 1.20,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'Cl': 1.75,
    'Br': 1.85,
    'P': 1.80,
    'S': 1.80,
    'As': 1.85,
    'B': 2.13,
    'Si': 2.10,
    'Se': 1.90,
}


atom_contrib = {
    table.GetAtomicNumber(s): 4. / 3. * pi * r ** 3
    for s, r in bondi_radii.items()
}


class VdwVolumeABC(Descriptor):
    r"""van der waals volume(ABC) descriptor.

    :returns: NaN when any atoms are non-compat_atoms

    References
        * :cite:`10.1021/jo034808o`
    """

    compat_atoms = tuple(bondi_radii)

    __slots__ = ()

    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return 'Vabc'

    def as_key(self):
        return self.__class__, ()

    def dependencies(self):
        return {
            'Nb': BondCount(),
            'NRa': RingCount(None, False, False, True, None),
            'NRA': RingCount(None, False, False, False, None),
        }

    def calculate(self, mol, Nb, NRa, NRA):
        ac = sum(
            atom_contrib.get(a.GetAtomicNum(), np.nan)
            for a in mol.GetAtoms()
        )

        return ac - 5.92 * Nb - 14.7 * NRa - 3.8 * NRA

    rtype = float
