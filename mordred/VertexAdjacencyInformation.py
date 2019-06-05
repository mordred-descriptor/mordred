import numpy as np

from ._base import Descriptor

__all__ = ("VertexAdjacencyInformation",)


class VertexAdjacencyInformation(Descriptor):
    r"""vertex adjacency information descriptor.

    .. math::
        {\rm VAdjMat} = 1 + \log_2(m)

    where :math:`m` is number of heavy-heavy bonds.

    :returns: :math:`m = 0`
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "vertex adjacency information"

    @classmethod
    def preset(cls, version):
        yield cls()

    explicit_hydrogens = False

    def __str__(self):
        return "VAdjMat"

    def parameters(self):
        return ()

    def calculate(self):
        m = sum(
            1
            for b in self.mol.GetBonds()
            if b.GetBeginAtom().GetAtomicNum() != 1
            and b.GetEndAtom().GetAtomicNum() != 1
        )

        with self.rethrow_zerodiv():
            return 1 + np.log2(m)

    rtype = float
