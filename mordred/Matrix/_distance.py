from .._base import Descriptor
from .._common import DistanceMatrix as D
from ._matrix_attributes import methods, get_method


class DistanceMatrix(Descriptor):
    r'''
    distance matrix descriptor

    Parameters:
        method(str): matrix aggregate method

    Returns:
        float: result
    '''

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        return map(cls, methods)

    @property
    def descriptor_key(self):
        return self.make_key(self.method)

    @property
    def dependencies(self):
        return dict(
            result=self.method.make_key(
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
        return '{}_D'.format(self.method.__name__)

    def __init__(self, method='SpMax'):
        self.method = get_method(method)

    def calculate(self, mol, result):
        return result
