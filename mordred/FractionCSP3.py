from ._base import Descriptor
from rdkit.Chem import AllChem as Chem

__all__ = (
    "FractionCSP3",
)


class FractionCSP3(Descriptor):
    r""" the fraction of C atoms that are SP3 hybridized.
    """

    __slots__ = ("_type",)

    @classmethod
    def preset(cls):
        yield cls()

    def description(self):
        return "the fraction of C atoms that are SP3 hybridized"

    def __str__(self):
        return self._type

    def __init__(self, type="FCSP3"):
        self._type = type

    def parameters(self):
        return ()

    def calculate(self):
        return Chem.CalcFractionCSP3(self.mol)

    rtype = int
