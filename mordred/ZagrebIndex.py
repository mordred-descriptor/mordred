from ._base import Descriptor
from ._graph_matrix import Valence

__all__ = ("ZagrebIndex",)


class ZagrebIndex(Descriptor):
    r"""Zagreb index descriptor.

    .. math::

        {}^\lambda M_1 = \sum_{atoms} d_i^\lambda

        {}^\lambda M_2 = \sum_{edges} \left(d_i \cdot d_j \right)^\lambda

    where
    :math:`d_i` is degree of i-th atom

    :type version: int
    :param version: Zagreb index version. 1 or 2.

    :type variable: int
    :param variable: lambda value.

    :returns: NaN when valence of any atoms are 0
    """

    since = "1.0.0"
    __slots__ = ("_version", "_variable")
    explicit_hydrogens = False

    def description(self):
        if self._variable == 1:
            return "Zagreb index (version {})".format(self._version)
        elif self._variable == -1:
            return "modified Zagreb index (version {})".format(self._version)
        else:
            return "Zagreb like index (lambda = {}, version {})".format(
                self._variable, self._version
            )

    @classmethod
    def preset(cls, version):
        return (cls(v, x) for x in [1, -1] for v in [1, 2])

    def __str__(self):
        if self._variable in {1, -1}:
            m = "" if self._variable == 1 else "m"
            return "{}Zagreb{}".format(m, self._version)

        return "Zagreb{}_{}".format(self._version, self._variable)

    def parameters(self):
        return self._version, self._variable

    def __init__(self, version=1, variable=1):
        assert version in {1, 2}
        self._version = version
        self._variable = variable

    def dependencies(self):
        return {"V": Valence(self.explicit_hydrogens)}

    def calculate(self, V):
        V = V.astype("float")

        if self._version == 1:
            with self.rethrow_zerodiv():
                return (V ** (self._variable * 2)).sum()
        else:
            return float(
                sum(
                    (V[b.GetBeginAtomIdx()] * V[b.GetEndAtomIdx()]) ** self._variable
                    for b in self.mol.GetBonds()
                )
            )

    rtype = float
