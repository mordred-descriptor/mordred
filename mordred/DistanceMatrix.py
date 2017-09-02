from ._base import Descriptor
from ._graph_matrix import DistanceMatrix as D
from ._matrix_attributes import methods, get_method

__all__ = ("DistanceMatrix",)


class DistanceMatrix(Descriptor):
    r"""distance matrix descriptor.

    :type type: str
    :param type: :ref:`matrix_aggregating_methods`
    """

    __slots__ = ("_type",)
    explicit_hydrogens = False

    def description(self):
        return "{} from distance matrix".format(self._type.description())

    @classmethod
    def preset(cls):
        return map(cls, methods)

    def __str__(self):
        return "{}_D".format(self._type.__name__)

    def parameters(self):
        return self._type,

    def __init__(self, type="SpMax"):
        self._type = get_method(type)

    def dependencies(self):
        return {
            "result": self._type(
                D(self.explicit_hydrogens),
                self.explicit_hydrogens,
                self.kekulize,
            ),
        }

    def calculate(self, result):
        return result

    rtype = float
