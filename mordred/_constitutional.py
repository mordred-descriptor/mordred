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

    :type prop: str or function
    :param prop: :ref:`atomic_properties`

    :rtype: float
    """

    require_connected = False

    @classmethod
    def preset(cls):
        return map(cls, _atomic_property.get_properties())

    _carbon = Chem.Atom(6)

    descriptor_keys = 'prop',

    def __init__(self, prop='v'):
        self.prop_name, self.prop = _atomic_property.getter(prop, self.explicit_hydrogens)
        if self.prop == _atomic_property.get_sanderson_en:
            self.prop_name = 'se'

    def __str__(self):
        return 'S{}'.format(self.prop_name)

    def calculate(self, mol):
        return sum(self.prop(a) / self.prop(self._carbon) for a in mol.GetAtoms())


class ConstitutionalMean(ConstitutionalSum):
    r"""mean of constitutional descriptor.

    .. math::
        M_p = \frac{S_p}{A}

    :type prop: str or function
    :param prop: :ref:`atomic_properties`

    :rtype: float
    """

    @classmethod
    def preset(cls):
        return map(cls, _atomic_property.get_properties())

    def __str__(self):
        return 'M{}'.format(self.prop_name)

    def dependencies(self):
        return dict(S=ConstitutionalSum(self.prop))

    def calculate(self, mol, S):
        return S / mol.GetNumAtoms()
