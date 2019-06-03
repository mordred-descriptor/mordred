from rdkit.Chem.Descriptors import MolWt, ExactMolWt

from ._base import Descriptor

__all__ = ("Weight",)


class Weight(Descriptor):
    r"""molecular weight descriptor.

    :type averaged: bool
    :param averaged: averaged by number of atom
    """

    def description(self):
        return "{}{}molecular weight".format(
            "averaged " if self._averaged else "", "exact " if self._exact else ""
        )

    since = "1.0.0"
    __slots__ = ("_averaged", "_exact")
    explicit_hydrogens = True

    @classmethod
    def preset(cls, version):
        yield cls(True, False)
        yield cls(True, True)

    def __str__(self):
        return "{}{}MW".format(
            "A" if self._averaged else "", "" if self._exact else "a"
        )

    def parameters(self):
        return self._exact, self._averaged

    def __init__(self, exact=True, averaged=False):
        self._averaged = averaged
        self._exact = exact

    def calculate(self):
        w = ExactMolWt(self.mol) if self._exact else MolWt(self.mol)
        if self._averaged:
            w /= self.mol.GetNumAtoms()

        return w

    rtype = float
