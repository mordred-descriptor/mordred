from .._base import Descriptor
from .._common import Eccentricity, Valence


class EccentricConnectivityIndex(Descriptor):
    explicit_hydrogens = False

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

    @property
    def descriptor_name(self):
        return 'ECIndex'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol, E, V):
        return (E.astype('int') * V).sum()
