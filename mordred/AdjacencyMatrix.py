from packaging.version import Version as StrictVersion

from ._base import Descriptor
from ._graph_matrix import AdjacencyMatrix as A
from ._matrix_attributes import SM1, methods, get_method

__all__ = ("AdjacencyMatrix",)

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version

_version_remove_SM1_A = StrictVersion("1.1.0")


class AdjacencyMatrix(Descriptor):
    r"""adjacency matrix descriptor.

    :type type: :py:class:`str`
    :param type: :ref:`matrix_aggregating_methods`
    """

    since = "1.0.0"
    __slots__ = ("_type",)
    explicit_hydrogens = False

    def description(self):
        return "{} of adjacency matrix".format(self._type.__name__)

    @classmethod
    def preset(cls, version):
        if version >= _version_remove_SM1_A:
            return (cls(m) for m in methods if m != SM1)
        else:
            return map(cls, methods)

    def __str__(self):
        return "{}_A".format(self._type.__name__)

    def parameters(self):
        return (self._type,)

    def __init__(self, type="SpMax"):
        self._type = get_method(type)

    def dependencies(self):
        return {
            "result": self._type(
                A(self.explicit_hydrogens), self.explicit_hydrogens, self.kekulize
            )
        }

    def calculate(self, result):
        return result

    rtype = float
