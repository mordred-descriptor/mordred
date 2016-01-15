from .._base import Descriptor
from .._common import Eccentricity, Valence


class EccentricConnectivityIndex(Descriptor):
    r'''
    eccentric connectivity index descriptor

    .. math::
        I_{\rm EC} = \sum^A_i {\boldsymbol E}{\boldsymbol V}

    where
    :math:`E` is eccentricity of atoms,
    :math:`V` is valences of atoms.

    Returns:
        int: eccentric connectivity index
    '''

    explicit_hydrogens = False
    descriptor_name = 'ECIndex'

    @property
    def dependencies(self):
        return dict(
            E=Eccentricity.make_key(
                self.explicit_hydrogens,
                False,
                False,
            ),
            V=Valence.make_key(
                self.explicit_hydrogens,
                False,
            ),
        )

    def calculate(self, mol, E, V):
        return (E.astype('int') * V).sum()


_descriptors = [EccentricConnectivityIndex]
__all__ = [d.__name__ for d in _descriptors]
