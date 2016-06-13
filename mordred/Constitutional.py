import numpy as np

from rdkit import Chem

from ._atomic_property import AtomicProperty, get_properties
from ._base import Descriptor


__all__ = ('ConstitutionalSum', 'ConstitutionalMean',)


class ConstitutionalSum(Descriptor):
    r"""sum of constitutional descriptor.

    .. math::
        S_p = \sum^A_{i=1} \frac{p_i}{p_{\rm C}}

    where
    :math:`p_i` is atomic property of i-th atom,
    :math:`p_{\rm C}` is atomic property of carbon

    :type prop: :py:class:`str` or :py:class:`function`
    :param prop: :ref:`atomic_properties`
    """
    __slots__ = ('_prop',)

    @classmethod
    def preset(cls):
        return map(cls, get_properties())

    _carbon = Chem.Atom(6)

    def as_key(self):
        return self.__class__, (self._prop,)

    def __init__(self, prop='v'):
        self._prop = AtomicProperty(self.explicit_hydrogens, prop)

    _prefix = 'S'

    def __str__(self):
        return '{}{}'.format(self._prefix, self._prop.as_argument)

    def dependencies(self):
        return {'P': self._prop}

    def calculate(self, P):
        C = self._prop.prop(self._carbon)
        return np.sum(P / C)

    rtype = float


class ConstitutionalMean(ConstitutionalSum):
    r"""mean of constitutional descriptor.

    .. math::
        M_p = \frac{S_p}{A}

    :type prop: :py:class:`str` or :py:class:`function`
    :param prop: :ref:`atomic_properties`

    :rtype: float
    """
    __slots__ = ('_prop',)
    _prefix = 'M'

    @classmethod
    def preset(cls):
        return map(cls, get_properties())

    def dependencies(self):
        return {'S': ConstitutionalSum(self._prop)}

    def calculate(self, S):
        return S / self.mol.GetNumAtoms()
