from ._base import Descriptor

import rdkit.Chem.rdMolDescriptors

__all__ = (
    "HeteroAtomsCount",
)


class HeteroAtomsCount(Descriptor):
    r"""hetero atom count descriptor.
    """

    __slots__ = ("_type",)

    @classmethod
    def preset(cls):
        yield cls()

    def description(self):
        return "number of hetero atoms"

    def __str__(self):
        return "n" + self._type

    def __init__(self, type="HAtoms"):
        self._type = type

    def parameters(self):
        return ()

    def calculate(self):
        return rdkit.Chem.rdMolDescriptors.CalcNumHeteroatoms(self.mol)

    rtype = int
