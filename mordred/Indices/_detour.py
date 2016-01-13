from .._base import Descriptor
from ..Matrix._detour import detour_matrix


class DetourIndex(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'DetourIndex'

    @property
    def dependencies(self):
        return dict(
            D=detour_matrix.make_key()
        )

    def calculate(self, mol, D):
        return int(0.5 * D.sum())
