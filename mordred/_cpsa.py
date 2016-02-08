import numpy as np

from ._base import Descriptor
from ._atomic_property import Rvdw, AtomicProperty
from .surface_area import SurfaceArea
from ._util import conformer_to_numpy, atoms_to_numpy


class CPSABase(Descriptor):
    require_3D = True

    @classmethod
    def preset(cls):
        yield cls()

    def __reduce_ex__(self, version):
        return self.__class__, ()

    def __str__(self):
        return self.__class__.__name__

    rtype = float


class VersionCPSABase(CPSABase):
    @classmethod
    def preset(cls):
        return map(cls, cls.versions)

    def __reduce_ex__(self, version):
        return self.__class__, (self._version,)

    versions = [1, 2, 3, 4, 5]

    def __init__(self, version=1):
        self._version = version

    def __str__(self):
        return '{}{}'.format(self.__class__.__name__, self._version)


class AtomicSurfaceArea(CPSABase):
    def __reduce_ex__(self, version):
        return self.__class__, (self._solvent_radius, self._level)

    def __init__(self, solvent_radius=1.4, level=5):
        self._solvent_radius = solvent_radius
        self._level = level

    def calculate(self, mol, conf):
        rs = atoms_to_numpy(lambda a: Rvdw[a.GetAtomicNum()] + self._solvent_radius, mol)

        ps = conformer_to_numpy(conf)

        sa = SurfaceArea(rs, ps, self._level)
        return np.array(sa.surface_area())


class TotalSurfaceArea(CPSABase):
    def dependencies(self):
        return {'ASA': AtomicSurfaceArea()}

    def calculate(self, mol, conf, ASA):
        return np.sum(ASA)

    rtype = float


class AtomicCharge(CPSABase):
    require_3D = False

    def dependencies(self):
        return {'charges': AtomicProperty(self.explicit_hydrogens, 'c')}

    def calculate(self, mol, charges):
        if not np.all(np.isfinite(charges)):
            return None

        else:
            return charges


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

    def calculate(self, mol, conf, SA, charges):
        if charges is None:
            return np.nan

        mask = self._mask(charges)

        if self._version == 1:
            f = 1.0
        elif self._version == 2:
            f = np.sum(charges[mask])
        elif self._version == 3:
            f = charges[mask]
        elif self._version == 4:
            f = np.sum(charges[mask]) / mol.GetNumAtoms()
        elif self._version == 5:
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
        return dict(
            PPSA=PPSA(self._version),
            PNSA=PNSA(self._version),
        )

    def calculate(self, mol, conf, PPSA, PNSA):
        return PPSA - PNSA


class FNSA(VersionCPSABase):
    r"""fractional charged partial negative surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    def dependencies(self):
        return dict(
            ASA=AtomicSurfaceArea(),
            SA=PNSA(self._version),
        )

    def calculate(self, mol, conf, SA, ASA):
        return SA / np.sum(ASA)


class FPSA(FNSA):
    r"""fractional charged partial positive surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    def dependencies(self):
        return dict(
            ASA=AtomicSurfaceArea(),
            SA=PPSA(self._version),
        )


class WNSA(FNSA):
    r"""surface weighted charged partial negative surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    def calculate(self, mol, conf, SA, ASA):
        return SA * np.sum(ASA) / 1000.0


class WPSA(FPSA):
    r"""surface weighted charged partial positive surface area descriptor.

    :type version: int
    :param version: one of :py:attr:`versions`
    """

    def calculate(self, mol, conf, SA, ASA):
        return SA * np.sum(ASA) / 1000.0


class RNCG(CPSABase):
    r"""relative negative charge descriptor."""

    require_3D = False

    @staticmethod
    def _mask(charges):
        return charges < 0.0

    def dependencies(self):
        return {'charges': AtomicCharge()}

    def calculate(self, mol, charges):
        if charges is None:
            return np.nan

        charges = charges[self._mask(charges)]
        Qmax = charges[np.argmax(np.abs(charges))]

        return Qmax / np.sum(charges)


class RPCG(RNCG):
    r"""relative positive charge descriptor."""

    @staticmethod
    def _mask(charges):
        return charges > 0.0


class RNCS(CPSABase):
    r"""relative negative charge surface area descriptor."""

    def dependencies(self):
        return dict(
            RCG=RNCG(),
            SA=AtomicSurfaceArea(),
            charges=AtomicCharge(),
        )

    @staticmethod
    def _mask(charges):
        return charges < 0

    def calculate(self, mol, conf, RCG, SA, charges):
        if charges is None:
            return np.nan

        mask = self._mask(charges)
        charges = charges[mask]

        SAmax = SA[mask][np.argmax(np.abs(charges))]

        return SAmax / RCG

    rtype = float


class RPCS(RNCS):
    r"""relative positive charge surface area descriptor."""

    @staticmethod
    def _mask(charges):
        return charges > 0

    def dependencies(self):
        return dict(
            RCG=RPCG(),
            SA=AtomicSurfaceArea(),
            charges=AtomicCharge(),
        )


class TASA(CPSABase):
    r"""total hydrophobic surface area descriptor."""

    @staticmethod
    def _mask(charges):
        return np.abs(charges) < 0.2

    def dependencies(self):
        return dict(
            SA=AtomicSurfaceArea(),
            charges=AtomicCharge(),
        )

    def calculate(self, mol, conf, SA, charges):
        if charges is None:
            return np.nan

        return np.sum(SA[self._mask(charges)])


class TPSA(TASA):
    r"""total polar surface area descriptor."""

    @staticmethod
    def _mask(charges):
        return np.abs(charges) >= 0.2


class RASA(CPSABase):
    r"""relative hydrophobic surface area descriptor."""

    def dependencies(self):
        return dict(
            TxSA=TASA(),
            SASA=AtomicSurfaceArea(),
        )

    def calculate(self, mol, conf, TxSA, SASA):
        return TxSA / np.sum(SASA)


class RPSA(RASA):
    r"""relative polar surface area descriptor."""

    def dependencies(self):
        return dict(
            TxSA=TPSA(),
            SASA=AtomicSurfaceArea(),
        )
