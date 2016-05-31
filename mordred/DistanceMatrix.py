from ._base import Descriptor
from ._graph_matrix import DistanceMatrix as D
from ._matrix_attributes import get_method, methods


__all__ = ('DistanceMatrix',)


class DistanceMatrix(Descriptor):
    r"""distance matrix descriptor.

    :type type: str
    :param type: :ref:`matrix_aggregating_methods`
    """

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        return map(cls, methods)

    def __str__(self):
        return '{}_D'.format(self._type.__name__)

    __slots__ = ('_type',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._type,)

    def __init__(self, type='SpMax'):
        self._type = get_method(type)

    def dependencies(self):
        return {
            'result': self._type(
                D(self.explicit_hydrogens),
                self.explicit_hydrogens,
                self.kekulize,
            )
        }

    def calculate(self, mol, result):
        return result

    rtype = float


if __name__ == '__main__':
    from .__main__ import submodule
    submodule()
