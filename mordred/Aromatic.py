from ._base import Descriptor

__all__ = ("AromaticAtomsCount", "AromaticBondsCount")


class AromaticBase(Descriptor):
    __slots__ = ()

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return self._name

    def parameters(self):
        return ()

    rtype = int


class AromaticAtomsCount(AromaticBase):
    r"""aromatic atoms count descriptor."""

    since = "1.0.0"
    __slots__ = ()
    _name = "nAromAtom"

    def description(self):
        return "aromatic atoms count"

    def calculate(self):
        return sum(1 for a in self.mol.GetAtoms() if a.GetIsAromatic())


class AromaticBondsCount(AromaticBase):
    r"""aromatic bonds count descriptor."""

    since = "1.0.0"
    __slots__ = ()
    _name = "nAromBond"

    def description(self):
        return "aromatic bonds count"

    def calculate(self):
        return sum(1 for b in self.mol.GetBonds() if b.GetIsAromatic())
