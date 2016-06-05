r"""charged partial surface area descriptor.

References
    * :cite:`10.1021/ac00220a013`
"""
import numpy as np

from ._base import Descriptor
from ._atomic_property import AtomicProperty, vdw_radii
from .surface_area import SurfaceArea
from ._util import atoms_to_numpy


__all__ = (
    'PNSA', 'PPSA',
    'DPSA',
    'FNSA', 'FPSA',
    'WNSA', 'WPSA',
    'RNCG', 'RPCG',
    'RNCS', 'RPCS',
    'TASA', 'TPSA',
    'RASA', 'RPSA',
)


class CPSABase(Descriptor):
    require_3D = True

    @classmethod
    def preset(cls):
        yield cls()

    def as_key(self):
        return self.__class__, ()

    def __str__(self):
        return self.__class__.__name__

    rtype = float


class VersionCPSABase(CPSABase):
    @classmethod
    def preset(cls):
        return map(cls, cls.versions)

    def as_key(self):
        return self.__class__, (self._version,)

    versions = [1, 2, 3, 4, 5]

    def __init__(self, version=1):
        self._version = version

    def __str__(self):
        return '{}{}'.format(self.__class__.__name__, self._version)


class AtomicSurfaceArea(CPSABase):
    def as_key(self):
        return self.__class__, (self._solvent_radius, self._level)

    def __init__(self, solvent_radius=1.4, level=5):
        self._solvent_radius = solvent_radius
        self._level = level

    def calculate(self):
        rs = atoms_to_numpy(lambda a: vdw_radii[a.GetAtomicNum()] + self._solvent_radius, self.mol)

        sa = SurfaceArea(rs, self.coord, self._level)
        return np.array(sa.surface_area())

    rtype = None


class TotalSurfaceArea(CPSABase):
    def dependencies(self):
        return {'ASA': AtomicSurfaceArea()}

    def calculate(self, ASA):
        return np.sum(ASA)


class AtomicCharge(CPSABase):
    require_3D = False

    def dependencies(self):
        return {'charges': AtomicProperty(self.explicit_hydrogens, 'c')}

    def calculate(self, charges):
        return charges

    rtype = None


class PNSA(VersionCPSABase):
    r"""partial negative surface area descriptor.

    .. math::

        {\rm PNSA}_1 = \sum_{a-} {\rm SA}_a^-

    where :math:`\sum_{a-}` means sum over negative charged atoms,
    :math:`{\rm SA}_a^-` is atomic partial surface area.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    def dependencies(self):
        return {
            'SA': AtomicSurfaceArea(),
            'charges': AtomicCharge(),
        }

    @staticmethod
    def _mask(charges):
        return charges < 0.0

    def calculate(self, SA, charges):
        mask = self._mask(charges)

        if self._version == 1:
            f = 1.0
        elif self._version == 2:
            f = np.sum(charges[mask])
        elif self._version == 3:
            f = charges[mask]
        elif self._version == 4:
            f = np.sum(charges[mask]) / self.mol.GetNumAtoms()
        elif self._version == 5:
            with self.rethrow_zerodiv():
                f = np.sum(charges[mask]) / np.sum(mask)

        return np.sum(f * SA[mask])


class PPSA(PNSA):
    r"""partial positive surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    @staticmethod
    def _mask(charges):
        return charges > 0.0


class DPSA(VersionCPSABase):
    r"""difference in charged partial surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    def dependencies(self):
        return {
            'PPSA': PPSA(self._version),
            'PNSA': PNSA(self._version),
        }

    def calculate(self, PPSA, PNSA):
        return PPSA - PNSA


class FNSA(VersionCPSABase):
    r"""fractional charged partial negative surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    def _SA(self):
        return PNSA(self._version)

    def dependencies(self):
        return {
            'ASA': AtomicSurfaceArea(),
            'SA': self._SA(),
        }

    def calculate(self, SA, ASA):
        return SA / np.sum(ASA)


class FPSA(FNSA):
    r"""fractional charged partial positive surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    def _SA(self):
        return PPSA(self._version)


class WxSAMixin(object):
    def calculate(self, SA, ASA):
        return SA * np.sum(ASA) / 1000.0


class WNSA(WxSAMixin, FNSA):
    r"""surface weighted charged partial negative surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    pass


class WPSA(WxSAMixin, FPSA):
    r"""surface weighted charged partial positive surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    pass


class RNCG(CPSABase):
    r"""relative negative charge descriptor."""

    require_3D = False

    @staticmethod
    def _mask(charges):
        return charges < 0.0

    def dependencies(self):
        return {'charges': AtomicCharge()}

    def calculate(self, charges):
        charges = charges[self._mask(charges)]
        if len(charges) == 0:
            return 0.0

        Qmax = charges[np.argmax(np.abs(charges))]

        return Qmax / np.sum(charges)


class RPCG(RNCG):
    r"""relative positive charge descriptor."""

    @staticmethod
    def _mask(charges):
        return charges > 0.0


class RNCS(CPSABase):
    r"""relative negative charge surface area descriptor."""

    _RCG = RNCG()

    def dependencies(self):
        return {
            'RCG': self._RCG,
            'SA': AtomicSurfaceArea(),
            'charges': AtomicCharge(),
        }

    @staticmethod
    def _mask(charges):
        return charges < 0

    def calculate(self, RCG, SA, charges):
        mask = self._mask(charges)
        charges = charges[mask]
        if len(charges) == 0:
            return 0.0

        SAmax = SA[mask][np.argmax(np.abs(charges))]

        return SAmax / RCG

    rtype = float


class RPCS(RNCS):
    r"""relative positive charge surface area descriptor."""

    @staticmethod
    def _mask(charges):
        return charges > 0

    _RCG = RPCG()


class TASA(CPSABase):
    r"""total hydrophobic surface area descriptor."""

    @staticmethod
    def _mask(charges):
        return np.abs(charges) < 0.2

    def dependencies(self):
        return {
            'SA': AtomicSurfaceArea(),
            'charges': AtomicCharge(),
        }

    def calculate(self, SA, charges):
        return np.sum(SA[self._mask(charges)])


class TPSA(TASA):
    r"""total polar surface area descriptor."""

    @staticmethod
    def _mask(charges):
        return np.abs(charges) >= 0.2


class RASA(CPSABase):
    r"""relative hydrophobic surface area descriptor."""

    _TxSA = TASA()

    def dependencies(self):
        return {
            'TxSA': self._TxSA,
            'SASA': AtomicSurfaceArea(),
        }

    def calculate(self, TxSA, SASA):
        return TxSA / np.sum(SASA)


class RPSA(RASA):
    r"""relative polar surface area descriptor."""

    _TxSA = TPSA()
