from rdkit.Chem.AllChem import CalcPBF

from ._base import Descriptor

__all__ = (
    "PBF",
)


class PBF(Descriptor):
    r"""PBF descriptor.
    """
    __slots__ = ()
    since = "1.1.0"

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
        return CalcPBF(self.mol)

    rtype = float
