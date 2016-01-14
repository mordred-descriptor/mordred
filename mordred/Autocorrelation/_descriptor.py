from .._base import Descriptor
from .. import _atomic_property
import numpy as np
from .._common import DistanceMatrix


class AutocorrelationBase(Descriptor):
    explicit_hydrogens = True

    prop = None

    @property
    def gasteiger_charges(self):
        return getattr(self.prop, 'gasteiger_charges', True)


class AVec(AutocorrelationBase):
    @property
    def descriptor_key(self):
        return self.make_key(self.prop)

    def __init__(self, prop):
        self.prop = prop

    def calculate(self, mol):
        return np.array([self.prop(a) for a in mol.GetAtoms()])


class CAVec(AutocorrelationBase):
    @property
    def descriptor_key(self):
        return self.make_key(self.prop)

    @property
    def dependencies(self):
        return dict(avec=AVec.make_key(self.prop))

    def __init__(self, prop):
        self.prop = prop

    def calculate(self, mol, avec):
        return avec - avec.mean()


class GMat(AutocorrelationBase):
    @property
    def descriptor_key(self):
        return self.make_key(self.distance)

    @property
    def dependencies(self):
        return dict(dmat=DistanceMatrix.make_key(
            self.explicit_hydrogens,
            False,
            False))

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
        elif self.distance == 0:
            return s
        else:
            return s / 2


class Autocorrelation(AutocorrelationBase):
    @property
    def descriptor_name(self):
        return '{}{}{}'.format(self.__class__.__name__, self.distance, self.prop_name)

    @property
    def descriptor_key(self):
        return self.__class__.make_key(self.prop, self.distance)

    def __init__(self, distance=0, prop='m'):
        if prop == 'c':
            self.prop_name = 'c'
            self.prop = _atomic_property.get_charge_explicitHs

        else:
            self.prop_name, self.prop = _atomic_property.getter(prop)

        self.distance = distance

    @property
    def _avec(self):
        return AVec.make_key(self.prop)

    @property
    def _cavec(self):
        return CAVec.make_key(self.prop)

    @property
    def _gmat(self):
        return GMat.make_key(self.distance)

    @property
    def _gsum(self):
        return GSum.make_key(self.distance)

    @property
    def _ATS(self):
        return ATS.make_key(self.distance, self.prop)

    @property
    def _ATSC(self):
        return ATSC.make_key(self.distance, self.prop)

    @property
    def _AATSC(self):
        return AATSC.make_key(self.distance, self.prop)

MAX_DISTANCE = 8


class ATS(Autocorrelation):
    r'''
    Autocorrelation of Topological Structure descriptor

    a.k.a. Moreau-Broto autocorrelation descriptor

    .. math::
        {\rm ATS}_0 = \sum^{A}_{i=1} {\boldsymbol w}_i^2

        {\rm ATS}_k = \frac{1}{2}
            \left(
                {\boldsymbol w}^{\rm T} \cdot
                {}^k{\boldsymbol B} \cdot
                {\boldsymbol w}
            \right)

        {}^k{\boldsymbol B} =
            \begin{cases}
                1 & (d_{ij} =    k) \\
                0 & (d_{ij} \neq k)
            \end{cases}

    where
    :math:`{\boldsymbol w}` is atomic property vector,
    :math:`d_{ij}` is graph distance.

    Parameters:
        distance(int): graph distance(:math:`k`)
        property(str, function): atomic property

    Returns:
        float: ATS value
    '''

    @classmethod
    def preset(cls):
        return (cls(d, a) for a in 'mvepis' for d in range(MAX_DISTANCE + 1))

    @property
    def dependencies(self):
        return dict(avec=self._avec, gmat=self._gmat)

    def calculate(self, mol, avec, gmat):
        if not gmat.any():
            return np.nan

        if self.distance == 0:
            return float((avec ** 2).sum())
        else:
            return 0.5 * avec.dot(gmat).dot(avec)


class AATS(ATS):
    r'''
    averaged ATS descriptor

    .. math::

        {\rm AATS}_k = \frac{{\rm ATS}_k}{\Delta_k}

    where
    :math:`\Delta_k` is number of vertex pairs at distance equal to :math:`k`.

    Parameters:
        parameters: see ATS

    Returns:
        float: AATS value
    '''

    @property
    def dependencies(self):
        return dict(ATS=self._ATS, gsum=self._gsum)

    def calculate(self, mol, ATS, gsum):
        return ATS / gsum


class ATSC(Autocorrelation):
    r'''
    centered ATS descriptor

    ATS with :math:`{\boldsymbol w}_{\rm c}` property

    .. math::
        {\boldsymbol w}_{\rm c} = {\boldsymbol w} - \bar{\boldsymbol w}

    Parameters:
        parameters: see ATS

    Returns:
        float: ATSC value

    '''

    @classmethod
    def preset(cls):
        return (cls(d, a) for a in 'cmvepis' for d in range(MAX_DISTANCE + 1))

    @property
    def dependencies(self):
        return dict(cavec=self._cavec, gmat=self._gmat)

    def calculate(self, mol, cavec, gmat):
        if not gmat.any():
            return np.nan

        if self.distance == 0:
            return float((cavec ** 2).sum())
        else:
            return 0.5 * cavec.dot(gmat).dot(cavec)


class AATSC(ATSC):
    r'''
    averaged ATSC descriptor

    .. math::

        {\rm AATSC}_k = \frac{{\rm ATSC}_k}{\Delta_k}

    where
    :math:`\Delta_k` is number of vertex pairs at distance equal to :math:`k`.

    Parameters:
        parameters: see ATS

    Returns:
        float: AATSC value
    '''

    @property
    def dependencies(self):
        return dict(ATSC=self._ATSC, gsum=self._gsum)

    def calculate(self, mol, ATSC, gsum):
        return ATSC / gsum


class MATS(Autocorrelation):
    r'''
    Moran coefficient descriptor

    .. math::

        {\rm MATS}_k = \frac{
            {\rm AATSC}_k
            }{
            \frac{1}{A} \cdot \sum {\boldsymbol w}_{\rm c}^2
            }

    Parameters:
        parameters: see ATS

    Returns:
        float: MATS value
    '''

    @classmethod
    def preset(cls):
        return (cls(d, a) for a in 'cmvepis' for d in range(1, MAX_DISTANCE + 1))

    @property
    def dependencies(self):
        return dict(avec=self._avec, AATSC=self._AATSC, cavec=self._cavec)

    def calculate(self, mol, avec, AATSC, cavec):
        return len(avec) * AATSC / (cavec ** 2).sum()


class GATS(MATS):
    r'''
    Geary coefficient descriptor

    Parameters:
        parameters: see ATS

    Returns:
        float: GATS value
    '''

    @property
    def dependencies(self):
        return dict(avec=self._avec, gmat=self._gmat, gsum=self._gsum, cavec=self._cavec)

    def calculate(self, mol, avec, gmat, gsum, cavec):
        if not gmat.any():
            return np.nan

        W = np.tile(avec, (len(avec), 1))
        n = (gmat * (W - W.T) ** 2).sum() / (4 * gsum)
        d = (cavec ** 2).sum() / (len(avec) - 1)
        return n / d

_descriptors = [ATS, AATS, ATSC, AATSC, MATS, GATS]
__all__ = [d.__name__ for d in _descriptors]
