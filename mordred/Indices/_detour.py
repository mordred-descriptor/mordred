from .._base import Descriptor
from ..Matrix._detour import detour_matrix


class DetourIndex(Descriptor):
    r'''
    detour index descriptor

    .. math::

        I_{\rm D} = \frac{1}{A} \sum^A_{i=1} \sum^A_{j=1} {\boldsymbol D}_{ij}

    where
    :math:`D` is detour matrix,
    :math:`A` is number of atoms.

    Returns:
        int: detour index
    '''

    explicit_hydrogens = False
    descriptor_name = 'DetourIndex'

    @property
    def dependencies(self):
        return dict(
            D=detour_matrix.make_key()
        )

    def calculate(self, mol, D):
        return int(0.5 * D.sum())
