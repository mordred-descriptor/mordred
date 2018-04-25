from __future__ import division

from ._base import Descriptor

__all__ = ("FragmentComplexity",)


class FragmentComplexity(Descriptor):
    r"""fragment complexity descriptor.

    .. math::
        {\rm fragCpx} = \left| B^2 - A^2 + A \right| + \frac{H}{100}

    where
    :math:`A` is number of atoms,
    :math:`B` is number of bonds,
    :math:`H` is number of hetero atoms

    References
        * :doi:`10.1021/ci050521b`

    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "fragment complexity"

    @classmethod
    def preset(cls, version):
        yield cls()

    explicit_hydrogens = False

    def parameters(self):
        return ()

    def __str__(self):
        return "fragCpx"

    def calculate(self):
        A = self.mol.GetNumAtoms()
        B = self.mol.GetNumBonds()
        H = sum(1 for a in self.mol.GetAtoms() if a.GetAtomicNum() != 6)
        return abs(B ** 2 - A ** 2 + A) + H / 100

    rtype = float
