from packaging.version import Version as StrictVersion

from rdkit.Chem import rdMolDescriptors

from ._base import Descriptor
from ._atomic_property import halogen

__all__ = ("AtomCount",)

_desc_conv = {
    "Atom": "all",
    "Bridgehead": "bridgehead",
    "HeavyAtom": "heavy",
    "Hetero": "hetero",
    "Spiro": "spiro",
    "X": "halogen",
}


_version_add_Nhetero = StrictVersion("1.1.0")


class AtomCount(Descriptor):
    r"""atom count descriptor.

    :type type: str
    :param type: type to count.

        * "Atom"
        * "HeavyAtom"
        * "Spiro"
        * "Bridgehead"
        * "Hetero"
        * "X" - all halogen
        * element symbol
    """

    since = "1.0.0"
    __slots__ = ("_type",)

    def description(self):
        return "number of {} atoms".format(_desc_conv.get(self._type, self._type))

    @classmethod
    def preset(cls, version):
        if version >= _version_add_Nhetero:
            t = [
                "Atom",
                "HeavyAtom",
                "Spiro",
                "Bridgehead",
                "Hetero",
                "H",
                "B",
                "C",
                "N",
                "O",
                "S",
                "P",
                "F",
                "Cl",
                "Br",
                "I",
                "X",
            ]
        else:
            t = [
                "Atom",
                "HeavyAtom",
                "Spiro",
                "Bridgehead",
                "H",
                "B",
                "C",
                "N",
                "O",
                "S",
                "P",
                "F",
                "Cl",
                "Br",
                "I",
                "X",
            ]

        return map(cls, t)

    @property
    def explicit_hydrogens(self):
        """Require explicit_hydrogens when type is "H" or "Atom"."""
        return self._type in {"H", "Atom"}

    def __str__(self):
        return "n" + self._type

    def parameters(self):
        return (self._type,)

    def __init__(self, type="Atom"):
        self._type = type

    def _calc_X(self):
        X = halogen
        return sum(a.GetAtomicNum() in X for a in self.mol.GetAtoms())

    def _calc(self):
        return sum(a.GetSymbol() == self._type for a in self.mol.GetAtoms())

    def _calc_all(self):
        return self.mol.GetNumAtoms()

    def calculate(self):
        if self._type == "X":
            return self._calc_X()
        elif self._type in ["Atom", "HeavyAtom"]:
            return self._calc_all()
        elif self._type == "Spiro":
            return rdMolDescriptors.CalcNumSpiroAtoms(self.mol)
        elif self._type == "Bridgehead":
            return rdMolDescriptors.CalcNumBridgeheadAtoms(self.mol)
        elif self._type == "Hetero":
            return rdMolDescriptors.CalcNumHeteroatoms(self.mol)
        else:
            return self._calc()

    rtype = int
