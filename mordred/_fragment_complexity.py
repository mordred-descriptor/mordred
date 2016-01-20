from ._base import Descriptor


class FragmentComplexity(Descriptor):
    r"""fragment complexity descriptor.

    .. math::
        {\rm fragCpx} = \left| B^2 - A^2 + A \right| + \frac{H}{100}

    where
    :math:`A` is number of atoms,
    :math:`B` is number of bonds,
    :math:`H` is number of hetero atoms

    :rtype: float

    References
        * :cite:`10.1021/ci050521b`
    """

    __slots__ = ()

    @classmethod
    def preset(cls):
        yield cls()

    explicit_hydrogens = False
    require_connected = False

    def __str__(self):
        return 'fragCpx'

    def calculate(self, mol):
        A = mol.GetNumAtoms()
        B = mol.GetNumBonds()
        H = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 6)
        return abs(B ** 2 - A ** 2 + A) + float(H) / 100
