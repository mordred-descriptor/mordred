import numpy as np

from rdkit import Chem

from ._atomic_property import AtomicProperty, get_properties, get_sanderson_en
from ._base import Descriptor


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

    @classmethod
    def preset(cls):
        return map(cls, get_properties())

    _carbon = Chem.Atom(6)

    __slots__ = ('_prop', '_prop_name',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._prop,)

    def __init__(self, prop='v'):
        self._prop = AtomicProperty(self.explicit_hydrogens, prop)

        if self._prop.prop == get_sanderson_en:
            self._prop_name = 'se'
        else:
            self._prop_name = str(self._prop)

    def __str__(self):
        return 'S{}'.format(self._prop_name)

    def dependencies(self):
        return {'P': self._prop}

    def calculate(self, mol, P):
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

    __slots__ = ('_prop_name', '_prop',)

    @classmethod
    def preset(cls):
        return map(cls, get_properties())

    def __str__(self):
        return 'M{}'.format(self._prop_name)

    def dependencies(self):
        return {'S': ConstitutionalSum(self._prop)}

    def calculate(self, mol, S):
        return S / mol.GetNumAtoms()
