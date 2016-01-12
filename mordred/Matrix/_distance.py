from .._base import Descriptor
from .._common import DistanceMatrix as D
from ._matrix_attributes import methods, get_method


class DistanceMatrix(Descriptor):
    explicit_hydrogens = False

    descriptor_defaults = [(m.__name__,) for m in methods]

    @property
    def descriptor_key(self):
        return self.make_key(self.method)

    @property
    def dependencies(self):
        return dict(
            result=get_method(self.method).make_key(
                D.make_key(
                    self.explicit_hydrogens,
                    False,
                    False,
                ),
                self.explicit_hydrogens,
                self.gasteiger_charges,
                self.kekulize,
            )
        )

    @property
    def descriptor_name(self):
        return '{}_D'.format(self.method)

    def __init__(self, method):
        self.method = method

    def calculate(self, mol, result):
        return result
