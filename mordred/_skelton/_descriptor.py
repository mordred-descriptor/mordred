from .._base import Descriptor


class NAME(Descriptor):
    explicit_hydrogens = False
    gasteiger_charges = True
    kekulize = True

    descriptor_defaults = []

    @property
    def dependencies(self):
        return dict()

    @property
    def descriptor_name(self):
        pass

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol):
        pass

_descriptors = []
__all__ = [d.__name__ for d in _descriptors]
