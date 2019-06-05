from rdkit.Chem.rdMolDescriptors import CalcPBF

from ._base import Descriptor

__all__ = ("PBF",)


class PBF(Descriptor):
    r"""PBF descriptor."""

    __slots__ = ()
    since = "1.1.2"
    require_3D = True

    @classmethod
    def preset(cls, version):
        yield cls()

    def description(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__

    def parameters(self):
        return ()

    def calculate(self):
        return CalcPBF(self.get_3D_mol())

    rtype = float
