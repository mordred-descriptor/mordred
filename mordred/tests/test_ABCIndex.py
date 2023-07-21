from rdkit import Chem
from numpy.testing import assert_almost_equal

from mordred.ABCIndex import ABCIndex, ABCGGIndex

# doi:10.2298/JSC150901093F
data = """
CC(C)CCCCCCC   6.58 6.49
CCC(C)CCCCCC   6.47 6.58
CC(C)(C)CCCCCC 6.84 6.82
CCC(C)(C)CCCCC 6.68 6.95
"""[
    1:-1
]


def test_ABC():
    abc = ABCIndex()
    abcgg = ABCGGIndex()

    for line in data.split("\n"):
        smi, dABC, dABCGG = line.strip().split()

        mol = Chem.MolFromSmiles(smi)
        assert assert_almost_equal(abc(mol), float(dABC), 2) is None
        assert assert_almost_equal(abcgg(mol), float(dABCGG), 2) is None
