from .._base import Descriptor
from .. import _atomic_property
from rdkit import Chem


class Sum(Descriptor):
    r'''
    sum of constitutional descriptor

    .. math::
        S_p = \sum^A_{i=1} \frac{p_i}{p_{\rm C}}

    where
    :math:`p_i` is atomic property of i-th atom,
    :math:`p_{\rm C}` is atomic property of carbon

    Parameters:
        prop(str, function): atomic property

    Returns:
        float: sum of constitutional
    '''

    @classmethod
    def preset(cls):
        return map(cls, ['v', 'se', 'pe', 'are', 'p', 'i'])

    _carbon = Chem.Atom(6)

    def __init__(self, prop='v'):
        self.prop_name, self.prop = _atomic_property.getter(prop)
        if self.prop == _atomic_property.get_sanderson_en:
            self.prop_name = 'se'

    @property
    def descriptor_name(self):
        return 'S{}'.format(self.prop_name)

    @property
    def descriptor_key(self):
        return self.make_key(self.prop)

    def calculate(self, mol):
        return sum(self.prop(a) / self.prop(self._carbon) for a in mol.GetAtoms())


class Mean(Sum):
    r'''
    mean of constitutional descriptor

    .. math::
        M_p = \frac{S_p}{A}

    Parameters:
        prop(str, function): atomic property

    Returns:
        float: mean of constitutional
    '''

    @classmethod
    def preset(cls):
        return map(cls, ['v', 'se', 'pe', 'are', 'p', 'i'])

    @property
    def descriptor_name(self):
        return 'M{}'.format(self.prop_name)

    @property
    def dependencies(self):
        return dict(S=Sum.make_key(self.prop))

    def calculate(self, mol, S):
        return S / mol.GetNumAtoms()


_descriptors = [Sum, Mean]
__all__ = [d.__name__ for d in _descriptors]
