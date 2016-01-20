from ._base import Descriptor
from ._detour_matrix import DetourMatrixCache


class DetourIndex(Descriptor):
    r"""detour index descriptor.

    .. math::

        I_{\rm D} = \frac{1}{A} \sum^A_{i=1} \sum^A_{j=1} {\boldsymbol D}_{ij}

    where
    :math:`D` is detour matrix,
    :math:`A` is number of atoms.

    :rtype: int
    """

    @classmethod
    def preset(cls):
        yield cls()

    explicit_hydrogens = False

    def __str__(self):
        return 'DetourIndex'

    def dependencies(self):
        return dict(
            D=DetourMatrixCache()
        )

    def calculate(self, mol, D):
        return int(0.5 * D.sum())
