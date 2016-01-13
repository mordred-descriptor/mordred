from .._base import Descriptor
from rdkit import Chem
from .. import _atomic_property
import numpy as np
from .._common import DistanceMatrix


class AutocorrelationBase(Descriptor):
    explicit_hydrogens = True

    attribute = None

    @property
    def gasteiger_charges(self):
        return getattr(self.attribute, 'gasteiger_charges', True)


class AVec(AutocorrelationBase):
    @property
    def descriptor_key(self):
        return self.make_key(self.attribute)

    def __init__(self, attribute):
        self.attribute = attribute

    def calculate(self, mol):
        return np.array([self.attribute(a) for a in mol.GetAtoms()])


class CAVec(AutocorrelationBase):
    @property
    def descriptor_key(self):
        return self.make_key(self.attribute)

    @property
    def dependencies(self):
        return dict(avec=AVec.make_key(self.attribute))

    def __init__(self, attribute):
        self.attribute = attribute

    def calculate(self, mol, avec):
        return avec - avec.mean()


class GMat(AutocorrelationBase):
    @property
    def descriptor_key(self):
        return self.make_key(self.distance)

    @property
    def dependencies(self):
        return dict(dmat=DistanceMatrix.make_key(True, False, False))

    def __init__(self, distance):
        self.distance = distance

    def calculate(self, mol, dmat):
        return dmat == self.distance


class GSum(AutocorrelationBase):
    @property
    def descriptor_key(self):
        return self.make_key(self.distance)

    @property
    def dependencies(self):
        return dict(gmat=GMat.make_key(self.distance))

    def __init__(self, distance):
        self.distance = distance

    def calculate(self, mol, gmat):
        s = gmat.sum()
        if s == 0:
            return np.nan
        else:
            return s


class Autocorrelation(AutocorrelationBase):
    @property
    def descriptor_name(self):
        return '{}{}{}'.format(self.__class__.__name__, self.distance, self.attr_name)

    @property
    def descriptor_key(self):
        return self.__class__.make_key(self.attribute, self.distance)

    def __init__(self, distance=0, attribute='m'):
        if attribute == 'c':
            self.attr_name = 'c'
            self.attribute = _atomic_property.get_charge_explicitHs

        else:
            self.attr_name, self.attribute = _atomic_property.getter(attribute)

        self.distance = distance

    @property
    def _avec(self):
        return AVec.make_key(self.attribute)

    @property
    def _cavec(self):
        return CAVec.make_key(self.attribute)

    @property
    def _gmat(self):
        return GMat.make_key(self.distance)

    @property
    def _gsum(self):
        return GSum.make_key(self.distance)

    @property
    def _ATS(self):
        return ATS.make_key(self.distance, self.attribute)

    @property
    def _ATSC(self):
        return ATSC.make_key(self.distance, self.attribute)

    @property
    def _AATSC(self):
        return AATSC.make_key(self.distance, self.attribute)

MAX_DISTANCE = 8


class ATS(Autocorrelation):
    @classmethod
    def preset(cls):
        return (cls(d, a) for a in 'mvepis' for d in range(MAX_DISTANCE + 1))

    @property
    def dependencies(self):
        return dict(avec=self._avec, gmat=self._gmat)

    def calculate(self, mol, avec, gmat):
        if not gmat.any():
            return np.nan

        r = float(avec.dot(gmat).dot(avec))
        return r if self.distance == 0 else r / 2.0


class AATS(ATS):
    @property
    def dependencies(self):
        return dict(ATS=self._ATS, gsum=self._gsum)

    def calculate(self, mol, ATS, gsum):
        r = ATS / gsum
        return r if self.distance == 0 else r * 2.0


class ATSC(Autocorrelation):
    @classmethod
    def preset(cls):
        return (cls(d, a) for a in 'cmvepis' for d in range(MAX_DISTANCE + 1))

    @property
    def dependencies(self):
        return dict(cavec=self._cavec, gmat=self._gmat)

    def calculate(self, mol, cavec, gmat):
        if not gmat.any():
            return np.nan

        r = float(cavec.dot(gmat).dot(cavec))
        return r if self.distance == 0 else r / 2.0


class AATSC(ATSC):
    @property
    def dependencies(self):
        return dict(ATSC=self._ATSC, gsum=self._gsum)

    def calculate(self, mol, ATSC, gsum):
        r = ATSC / gsum
        return r if self.distance == 0 else r * 2.0


class MATS(Autocorrelation):
    @classmethod
    def preset(cls):
        return (cls(d, a) for a in 'cmvepis' for d in range(1, MAX_DISTANCE + 1))

    @property
    def dependencies(self):
        return dict(avec=self._avec, AATSC=self._AATSC, cavec=self._cavec)

    def calculate(self, mol, avec, AATSC, cavec):
        return len(avec) * AATSC / (cavec ** 2).sum()


class GATS(MATS):
    @property
    def dependencies(self):
        return dict(avec=self._avec, gmat=self._gmat, gsum=self._gsum, cavec=self._cavec)

    def calculate(self, mol, avec, gmat, gsum, cavec):
        if not gmat.any():
            return np.nan

        W = np.tile(avec, (len(avec), 1))
        n = (gmat * (W - W.T) ** 2).sum() / (2 * gsum)
        d = (cavec ** 2).sum() / (len(avec) - 1)
        return n / d

_descriptors = [ATS, AATS, ATSC, AATSC, MATS, GATS]
__all__ = [d.__name__ for d in _descriptors]
