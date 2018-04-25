from ._base import Descriptor
from ._graph_matrix import Valence, Eccentricity

__all__ = ("EccentricConnectivityIndex",)


class EccentricConnectivityIndex(Descriptor):
    r"""eccentric connectivity index descriptor.

    .. math::
        I_{\rm EC} = \sum^A_i {\boldsymbol E}{\boldsymbol V}

    where
    :math:`E` is eccentricity of atoms,
    :math:`V` is valences of atoms.
    """

    since = "1.0.0"
    __slots__ = ()
    explicit_hydrogens = False

    def description(self):
        return "eccentric connectivity index"

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return "ECIndex"

    def parameters(self):
        return ()

    def dependencies(self):
        return {
            "E": Eccentricity(self.explicit_hydrogens),
            "V": Valence(self.explicit_hydrogens),
        }

    def calculate(self, E, V):
        return int((E.astype("int") * V).sum())

    rtype = int
