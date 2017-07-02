from six import string_types, integer_types
from numpy import product

from ._base import Descriptor
from ._graph_matrix import Valence, DistanceMatrix
from ._atomic_property import table

__all__ = ("MolecularDistanceEdge",)


class MolecularDistanceEdge(Descriptor):
    r"""molecular distance edge descriptor.

    :type valence1: int
    :param valence1: valence of first atom

    :type valence2: int
    :param valence2: valence of second atom

    :type element: str or int
    :param element: atomic symbol or atomic number

    :returns: NaN when :math:`N_{\rm MDE} = 0`
    """

    __slots__ = ("_valence1", "_valence2", "_atomic_num",)
    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        return (
            cls(a, b, e)
            for e in [6, 8, 7]
            for a in range(1, 11 - e)
            for b in range(a, 11 - e)
        )

    def __str__(self):
        return "MDE{}-{}{}".format(
            table.GetElementSymbol(self._atomic_num),
            self._valence1,
            self._valence2,
        )

    def parameters(self):
        return self._valence1, self._valence2, self._atomic_num

    def __init__(self, valence1=1, valence2=1, element="C"):
        self._valence1 = min(valence1, valence2)
        self._valence2 = max(valence1, valence2)
        if isinstance(element, integer_types):
            self._atomic_num = element
        elif isinstance(element, string_types):
            self._atomic_num = table.GetAtomicNumber(element)
        else:
            raise ValueError("element must be atomic number or atomic symbol")

    def dependencies(self):
        return {
            "D": DistanceMatrix(self.explicit_hydrogens),
            "V": Valence(self.explicit_hydrogens),
        }

    def calculate(self, D, V):
        N = self.mol.GetNumAtoms()
        Dv = [
            D[i, j]
            for i in range(N)
            for j in range(i + 1, N)
            if (V[i] == self._valence1 and V[j] == self._valence2) or
            (V[j] == self._valence1 and V[i] == self._valence2)
            if self.mol.GetAtomWithIdx(i).GetAtomicNum() ==
            self.mol.GetAtomWithIdx(j).GetAtomicNum() ==
            self._atomic_num
        ]
        n = len(Dv)

        with self.rethrow_zerodiv():
            dx = product(Dv) ** (1. / (2. * n))

        return n / (dx ** 2)

    rtype = float
