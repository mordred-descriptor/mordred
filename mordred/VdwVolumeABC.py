from math import pi

from ._base import Descriptor
from .BondCount import BondCount
from .RingCount import RingCount
from ._atomic_property import GetAtomicNumber

__all__ = ("VdwVolumeABC",)


bondi_radii = {
    "H": 1.20,
    "C": 1.70,  # noqa: S001
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,  # noqa: S001
    "Cl": 1.75,  # noqa: S001
    "Br": 1.85,  # noqa: S001
    "P": 1.80,
    "S": 1.80,
    "As": 1.85,  # noqa: S001
    "B": 2.13,
    "Si": 2.10,
    "Se": 1.90,  # noqa: S001
}


atom_contrib = {
    GetAtomicNumber(s): 4.0 / 3.0 * pi * r ** 3 for s, r in bondi_radii.items()
}


class VdwVolumeABC(Descriptor):
    r"""van der waals volume(ABC) descriptor.

    :returns: NaN when any atoms are non-compat_atoms

    References
        * :doi:`10.1021/jo034808o`

    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "ABC van der waals volume"

    compat_atoms = tuple(bondi_radii)

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return "Vabc"

    def parameters(self):
        return ()

    def dependencies(self):
        return {
            "NRA": RingCount(None, False, False, False, None),
            "NRa": RingCount(None, False, False, True, None),
            "Nb": BondCount(),
        }

    def calculate(self, Nb, NRa, NRA):
        try:
            ac = sum(atom_contrib[a.GetAtomicNum()] for a in self.mol.GetAtoms())
        except KeyError:
            self.fail(ValueError("unknown atom type"))

        return ac - 5.92 * Nb - 14.7 * NRa - 3.8 * NRA

    rtype = float
