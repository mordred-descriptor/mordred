from rdkit import Chem
from numpy.testing import assert_almost_equal

from mordred.LogS import LogS


def test_LogS1():
    mol = Chem.MolFromSmiles("Oc1ccccc1OC")
    yield assert_almost_equal, LogS()(mol), -1.91, 2


def test_LogS2():
    mol = Chem.MolFromSmiles("c1ccccc1c2ccccn2")
    yield assert_almost_equal, LogS()(mol), -4.53, 2
