# Based on the code from https://github.com/rdkit/rdkit/blob/master/Contrib/PBF/pbf.py

from ._base import Descriptor

__all__ = (
    "SP3A2C",
)


class SP3A2C(Descriptor):
    r"""SP3A2C descriptor.
    """

    @classmethod
    def preset(cls):
        yield cls()

    def description(self):
        return self.__class__.__name__

    def __str__(self):
        return self.__class__.__name__

    def __init__(self, type="PBF"):
        self._type = type

    def parameters(self):
        return ()

    def calculate(self):
        acyclic = 0
        cyclic = 0
        for x in self.mol.GetAtoms():
            if x.GetSymbol() == 'C' and str(x.GetHybridization()) == 'SP3':
                if x.IsInRing():
                    cyclic += 1
                else:
                    acyclic += 1

        sp3 = cyclic / (float(cyclic + acyclic))
        return sp3

    rtype = float
