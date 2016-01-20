import numpy as np

from . import _atomic_property
from ._base import Descriptor
from ._common import DistanceMatrix


class AutocorrelationBase(Descriptor):
    explicit_hydrogens = True
    prop = None

    @property
    def gasteiger_charges(self):
        return getattr(self.prop, 'gasteiger_charges', False)

    @property
    def require_connected(self):
        return getattr(self.prop, 'require_connected', False)

    def __str__(self):
        return '{}{}{}'.format(
            self.__class__.__name__,
            self.order,
            self.prop_name
        )

    descriptor_keys = 'order', 'prop'

    def __init__(self, order=0, prop='m'):
        self.prop_name, self.prop = _atomic_property.getter(prop, self.explicit_hydrogens)
        self.order = order

    @property
    def _avec(self):
        return AVec(self.prop)

    @property
    def _cavec(self):
        return CAVec(self.prop)

    @property
    def _gmat(self):
        return GMat(self.order)

    @property
    def _gsum(self):
        return GSum(self.order)

    @property
    def _ATS(self):
        return ATS(self.order, self.prop)

    @property
    def _ATSC(self):
        return ATSC(self.order, self.prop)

    @property
    def _AATSC(self):
        return AATSC(self.order, self.prop)


class AVec(AutocorrelationBase):
    descriptor_keys = 'prop',

    def __init__(self, prop):
        self.prop = prop

    def calculate(self, mol):
        return np.array([self.prop(a) for a in mol.GetAtoms()])


class CAVec(AutocorrelationBase):
    descriptor_keys = 'prop',

    def __init__(self, prop):
        self.prop = prop

    def dependencies(self):
        return dict(avec=AVec(self.prop))

    def calculate(self, mol, avec):
        return avec - avec.mean()


class GMat(AutocorrelationBase):
    descriptor_keys = 'order',

    def __init__(self, order):
        self.order = order

    def dependencies(self):
        return dict(
            dmat=DistanceMatrix(
                self.explicit_hydrogens,
                False,
                False,
            )
        )

    def calculate(self, mol, dmat):
        return dmat == self.order


class GSum(AutocorrelationBase):
    descriptor_keys = 'order',

    def __init__(self, order):
        self.order = order

    def dependencies(self):
        return dict(gmat=GMat(self.order))

    def calculate(self, mol, gmat):
        s = gmat.sum()

        if self.order == 0:
            return s
        else:
            return s / 2


MAX_DISTANCE = 8


class ATS(AutocorrelationBase):
    r"""Autocorrelation of Topological Structure descriptor.

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
    :math:`d_{ij}` is graph distance(smallest number of bonds between atom i and j).

    :type order: int
    :param order: order(:math:`k`)

    :type property: str, function
    :param property: :ref:`atomic_properties`

    :type: float
    """

    @classmethod
    def preset(cls):
        return (cls(d, a)
                for a in _atomic_property.get_properties(istate=True)
                for d in range(MAX_DISTANCE + 1))

    def dependencies(self):
        return dict(avec=self._avec, gmat=self._gmat)

    def calculate(self, mol, avec, gmat):
        if self.order == 0:
            return float((avec ** 2).sum())

        return 0.5 * avec.dot(gmat).dot(avec)


class AATS(ATS):
    r"""averaged ATS descriptor.

    .. math::

        {\rm AATS}_k = \frac{{\rm ATS}_k}{\Delta_k}

    where
    :math:`\Delta_k` is number of vertex pairs at order equal to :math:`k`.

    :Parameters: see :py:class:`ATS`

    :rtype: float
    """

    def dependencies(self):
        return dict(ATS=self._ATS, gsum=self._gsum)

    def calculate(self, mol, ATS, gsum):
        return ATS / (gsum or np.nan)


class ATSC(AutocorrelationBase):
    r"""centered ATS descriptor.

    ATS with :math:`{\boldsymbol w}_{\rm c}` property

    .. math::
        {\boldsymbol w}_{\rm c} = {\boldsymbol w} - \bar{\boldsymbol w}

    :Parameters: see ATS

    :rtype: float
    """

    @classmethod
    def preset(cls):
        return (cls(d, a)
                for a in _atomic_property.get_properties(charge=True, istate=True)
                for d in range(MAX_DISTANCE + 1))

    def dependencies(self):
        return dict(cavec=self._cavec, gmat=self._gmat)

    def calculate(self, mol, cavec, gmat):
        if self.order == 0:
            return float((cavec ** 2).sum())

        return 0.5 * cavec.dot(gmat).dot(cavec)


class AATSC(ATSC):
    r"""averaged ATSC descriptor.

    .. math::

        {\rm AATSC}_k = \frac{{\rm ATSC}_k}{\Delta_k}

    where
    :math:`\Delta_k` is number of vertex pairs at order equal to :math:`k`.

    :Parameters: see :py:class:`ATS`

    :rtype: float
    """

    def dependencies(self):
        return dict(ATSC=self._ATSC, gsum=self._gsum)

    def calculate(self, mol, ATSC, gsum):
        return ATSC / (gsum or np.nan)


class MATS(AutocorrelationBase):
    r"""Moran coefficient descriptor.

    .. math::

        {\rm MATS}_k = \frac{
            {\rm AATSC}_k
            }{
            \frac{1}{A} \cdot \sum {\boldsymbol w}_{\rm c}^2
            }

    :Parameters: see :py:class:`ATS`

    :rtype: float
    """

    @classmethod
    def preset(cls):
        return (cls(d, a)
                for a in _atomic_property.get_properties(charge=True, istate=True)
                for d in range(1, MAX_DISTANCE + 1))

    def dependencies(self):
        return dict(avec=self._avec, AATSC=self._AATSC, cavec=self._cavec)

    def calculate(self, mol, avec, AATSC, cavec):
        return len(avec) * AATSC / (cavec ** 2).sum()


class GATS(MATS):
    r"""Geary coefficient descriptor.

    :Parameters: see :py:class:`ATS`

    :rtype: float
    """

    def dependencies(self):
        return dict(avec=self._avec, gmat=self._gmat, gsum=self._gsum, cavec=self._cavec)

    def calculate(self, mol, avec, gmat, gsum, cavec):
        W = np.tile(avec, (len(avec), 1))
        n = (gmat * (W - W.T) ** 2).sum() / (4 * (gsum or np.nan))
        d = (cavec ** 2).sum() / (len(avec) - 1)
        return n / d
