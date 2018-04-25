r"""Wildman-Crippen LogP/MR descriptor.

References
    * :doi:`10.1021/ci990307l`

"""
from rdkit.Chem import Crippen

from ._base import Descriptor

__all__ = ("SLogP", "SMR")


class WildmanCrippenBase(Descriptor):
    __slots__ = ()

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return self.__class__.__name__

    def parameters(self):
        return ()

    explicit_hydrogens = False

    rtype = float


class SLogP(WildmanCrippenBase):
    r"""Wildman-Crippen LogP descriptor(rdkit wrapper)."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "Wildman-Crippen LogP"

    def calculate(self):
        return Crippen.MolLogP(self.mol)


class SMR(WildmanCrippenBase):
    r"""Wildman-Crippen MR descriptor(rdkit wrapper)."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "Wildman-Crippen MR"

    def calculate(self):
        return Crippen.MolMR(self.mol)
