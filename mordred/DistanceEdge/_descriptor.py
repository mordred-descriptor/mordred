from .._base import Descriptor
from .._common import DistanceMatrix, Valence
from six import integer_types, string_types
from rdkit import Chem
from numpy import product, nan


table = Chem.GetPeriodicTable()


class DistanceEdge(Descriptor):
    explicit_hydrogens = False

    descriptor_defaults = [
        (1, 1, 6), (1, 2, 6), (1, 3, 6), (1, 4, 6),
        (2, 2, 6), (2, 3, 6), (2, 4, 6),
        (3, 3, 6), (3, 4, 6),
        (4, 4, 6),
        (1, 1, 8), (1, 2, 8),
        (2, 2, 8),
        (1, 1, 7), (1, 2, 7), (1, 3, 7),
        (2, 2, 7), (2, 3, 7),
        (3, 3, 7),
    ]

    @property
    def dependencies(self):
        return dict(
            D=DistanceMatrix.make_key(
                self.explicit_hydrogens,
                False,
                False,
            ),
            V=Valence.make_key(
                False,
                False,
            )
        )

    def __init__(self, valence1=1, valence2=1, element='C'):
        self.valence1 = min(valence1, valence2)
        self.valence2 = max(valence1, valence2)
        if isinstance(element, integer_types):
            self.atomic_num = element
        elif isinstance(element, string_types):
            self.atomic_num = table.GetAtomicNumber(element)
        else:
            raise ValueError('element must be atomic number or atomic symbol')

    @property
    def descriptor_name(self):
        return 'MDE{}-{}{}'.format(
            table.GetElementSymbol(self.atomic_num),
            self.valence1,
            self.valence2,
        )

    @property
    def descriptor_key(self):
        return self.make_key(
            self.valence1, self.valence2, self.atomic_num
        )

    def calculate(self, mol, D, V):
        N = mol.GetNumAtoms()
        Dv = [
            D[i, j]
            for i in range(N)
            for j in range(i + 1, N)
            if (V[i] == self.valence1 and V[j] == self.valence2)
            or (V[j] == self.valence1 and V[i] == self.valence2)
            if mol.GetAtomWithIdx(i).GetAtomicNum() ==
            mol.GetAtomWithIdx(j).GetAtomicNum() ==
            self.atomic_num
        ]
        n = len(Dv)
        if n == 0:
            return nan

        dx = product(Dv) ** (1. / (2. * n))

        return n / (dx ** 2)

_descriptors = [DistanceEdge]
__all__ = [d.__name__ for d in _descriptors]
