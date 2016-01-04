from .._base import Descriptor
from .. import _atomic_property
from rdkit import Chem


class Sum(Descriptor):
    descriptor_defaults = [('v',), ('se',), ('pe',), ('are',), ('p',), ('i',)]

    _carbon = Chem.Atom(6)

    def __init__(self, attribute='v'):
        self.attr_name, self.attribute = _atomic_property.getter(attribute)
        if self.attribute == _atomic_property.get_sanderson_en:
            self.attr_name = 'se'

    @property
    def descriptor_name(self):
        return 'S{}'.format(self.attr_name)

    @property
    def descriptor_key(self):
        return self.make_key(self.attribute)

    def calculate(self, mol):
        return sum((self.attribute(a) / self.attribute(self._carbon) for a in mol.GetAtoms()))


class Mean(Sum):
    descriptor_defaults = [('v',), ('se',), ('pe',), ('are',), ('p',), ('i',)]

    @property
    def descriptor_name(self):
        return 'M{}'.format(self.attr_name)

    @property
    def dependencies(self):
        return dict(S=Sum.make_key(self.attribute))

    def calculate(self, mol, S):
        return S / mol.GetNumAtoms()


_descriptors = [Sum, Mean]
__all__ = [d.__name__ for d in _descriptors]
