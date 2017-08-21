from rdkit.Chem.Descriptors import ExactMolWt

from ._base import Descriptor

__all__ = (
    "Weight",
)


class Weight(Descriptor):
    r"""molecular weight descriptor.

    :type averaged: bool
    :param averaged: averaged by number of atom
    """

    def description(self):
        return "{}molecular weight".format("averaged " if self._averaged else "")

    __slots__ = ("_averaged",)
    explicit_hydrogens = True

    @classmethod
    def preset(cls):
        yield cls(False)
        yield cls(True)

    def __str__(self):
        return "AMW" if self._averaged else "MW"

    def parameters(self):
        return self._averaged,

    def __init__(self, averaged=False):
        self._averaged = averaged

    def calculate(self):
        w = ExactMolWt(self.mol)
        if self._averaged:
            w /= self.mol.GetNumAtoms()

        return w

    rtype = float
