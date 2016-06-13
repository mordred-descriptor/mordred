from math import pi

from ._base import Descriptor
from .BondCount import BondCount
from .RingCount import RingCount
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
    __slots__ = ()

    compat_atoms = tuple(bondi_radii)

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

    def calculate(self, Nb, NRa, NRA):
        try:
            ac = sum(
                atom_contrib[a.GetAtomicNum()]
                for a in self.mol.GetAtoms()
            )
        except KeyError:
            self.fail(ValueError('unknown atom type'))

        return ac - 5.92 * Nb - 14.7 * NRa - 3.8 * NRA

    rtype = float
