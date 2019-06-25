from rdkit.Chem.QED import qed

from ._base import Descriptor

__all__ = (
    "QED",
)


class QED(Descriptor):
    r"""QED descriptor."""

    __slots__ = ()
    since = "1.1.2"
    explicit_hydrogens = False

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
        return qed(self.mol)

    rtype = float
