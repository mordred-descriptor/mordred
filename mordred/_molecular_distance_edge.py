from ._base import Descriptor
from ._common import DistanceMatrix, Valence
from six import integer_types, string_types
from rdkit import Chem
from numpy import product, nan


table = Chem.GetPeriodicTable()


class MolecularDistanceEdge(Descriptor):
    '''
    molecular distance edge descriptor

    Parameters:
        valence1(int): valence of first atom
        valence2(int): valence of second atom

        element(str, int): atomic symbol or atomic number

    Returns:
        float: MDE value
    '''

    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        return (
            cls(a, b, e)
            for e in [6, 8, 7]
            for a in range(1, 11 - e)
            for b in range(a, 11 - e)
        )

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
