from ._base import Descriptor
from ._graph_matrix import Eccentricity, Valence


__all__ = ('EccentricConnectivityIndex',)


class EccentricConnectivityIndex(Descriptor):
    r"""eccentric connectivity index descriptor.

    .. math::
        I_{\rm EC} = \sum^A_i {\boldsymbol E}{\boldsymbol V}

    where
    :math:`E` is eccentricity of atoms,
    :math:`V` is valences of atoms.
    """

    __slots__ = ()

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return 'ECIndex'

    def __reduce_ex__(self, version):
        return self.__class__, ()

    def dependencies(self):
        return {
            'E': Eccentricity(self.explicit_hydrogens),
            'V': Valence(self.explicit_hydrogens),
        }

    def calculate(self, mol, E, V):
        return int((E.astype('int') * V).sum())

    rtype = int
