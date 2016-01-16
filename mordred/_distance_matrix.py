from ._base import Descriptor
from ._common import DistanceMatrix as D
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

    def __str__(self):
        return '{}_D'.format(self.method.__name__)

    descriptor_keys = 'method',

    def __init__(self, method='SpMax'):
        self.method = get_method(method)

    @property
    def dependencies(self):
        return dict(
            result=self.method(
                D(
                    self.explicit_hydrogens,
                    False,
                    False,
                ),
                self.explicit_hydrogens,
                self.gasteiger_charges,
                self.kekulize,
            )
        )

    def calculate(self, mol, result):
        return result
