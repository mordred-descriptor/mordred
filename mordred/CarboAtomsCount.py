from ._base import Descriptor

__all__ = (
    "CarbonAtomsCount",
)


class CarbonAtomsCount(Descriptor):
    r"""carbon atom count descriptor.
    """

    __slots__ = ("_type",)

    @classmethod
    def preset(cls):
        yield cls()

    def description(self):
        return "number of carbon atoms"

    def __str__(self):
        return "n" + self._type

    def __init__(self, type="CAtoms"):
        self._type = type

    def parameters(self):
        return ()

    def calculate(self):
        num_carbon = 0
        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                num_carbon += 1
        return num_carbon

    rtype = int
