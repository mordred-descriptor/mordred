from .._base import Descriptor

__all__ = []


class NAME(Descriptor):
    @property
    def descriptor_key(self):
        pass

    def calculate(self, mol):
        pass

_descriptors = []
