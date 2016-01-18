from ._base import Descriptor
from ._detour_matrix import detour_matrix


class DetourIndex(Descriptor):
    r'''
    detour index descriptor

    .. math::

        I_{\rm D} = \frac{1}{A} \sum^A_{i=1} \sum^A_{j=1} {\boldsymbol D}_{ij}

    where
    :math:`D` is detour matrix,
    :math:`A` is number of atoms.

    :rtype: int
    '''

    explicit_hydrogens = False

    def __str__(self):
        return 'DetourIndex'

    @property
    def dependencies(self):
        return dict(
            D=detour_matrix()
        )

    def calculate(self, mol, D):
        return int(0.5 * D.sum())
