from .._base import Descriptor
from .._common import Eccentricity, Valence


class EccentricConnectivityIndex(Descriptor):
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
