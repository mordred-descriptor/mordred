from ._base import Descriptor
from ._common import DistanceMatrix as D
from ._matrix_attributes import methods, get_method


class DistanceMatrix(Descriptor):
    r'''
    distance matrix descriptor

    Parameters:
        type(str): matrix aggregating method

    Returns:
        float: result
    '''

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        return map(cls, methods)

    def __str__(self):
        return '{}_D'.format(self.type.__name__)

    descriptor_keys = 'type',

    def __init__(self, type='SpMax'):
        self.type = get_method(type)

    @property
    def dependencies(self):
        return dict(
            result=self.type(
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
