from rdkit import Chem

from . import _atomic_property
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
        return map(cls, _atomic_property.get_properties())

    _carbon = Chem.Atom(6)

    __slots__ = ('_prop_name', '_prop',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._prop,)

    def __init__(self, prop='v'):
        self._prop_name, self._prop = _atomic_property.getter(prop, self.explicit_hydrogens)
        if self._prop == _atomic_property.get_sanderson_en:
            self._prop_name = 'se'

    def __str__(self):
        return 'S{}'.format(self._prop_name)

    def calculate(self, mol):
        return sum(float(self._prop(a)) / self._prop(self._carbon) for a in mol.GetAtoms())

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
        return map(cls, _atomic_property.get_properties())

    def __str__(self):
        return 'M{}'.format(self._prop_name)

    def dependencies(self):
        return dict(S=ConstitutionalSum(self._prop))

    def calculate(self, mol, S):
        return S / mol.GetNumAtoms()
