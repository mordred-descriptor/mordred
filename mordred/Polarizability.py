from ._base import Descriptor
from ._atomic_property import polarizability78, polarizability94

__all__ = ("APol", "BPol")


class PolarizabilityBase(Descriptor):
    __slots__ = ("_use78",)

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return self.__class__.__name__.lower() + ("78" if self._use78 else "")

    def parameters(self):
        return (self._use78,)

    def __init__(self, use78=False):
        self._use78 = use78

    def _get_table(self):
        return polarizability78 if self._use78 else polarizability94

    rtype = float


class APol(PolarizabilityBase):
    r"""atomic polarizability descriptor.

    :type use78: bool
    :param use78: use old atomic polarizability data
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "atomic polarizability"

    def calculate(self):
        table = self._get_table()
        return sum(table[a.GetAtomicNum()] for a in self.mol.GetAtoms())


class BPol(PolarizabilityBase):
    r"""bond polarizability descriptor.

    :type use78: bool
    :param use78: use old atomic polarizability data
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "bond polarizability"

    def calculate(self):
        table = self._get_table()

        def bond_pol(bond):
            a = bond.GetBeginAtom().GetAtomicNum()
            b = bond.GetEndAtom().GetAtomicNum()
            return abs(table[a] - table[b])

        return float(sum(bond_pol(b) for b in self.mol.GetBonds()))
