from .._base import Descriptor
from ..Matrix._detour import detour_matrix


class Detour(Descriptor):
    explicit_hydrogens = False

    @property
    def dependencies(self):
        return dict(
            D=detour_matrix.make_key()
        )

    @property
    def descriptor_name(self):
        return 'Dt'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol, D):
        return int(0.5 * D.sum())
