from mordred.WildmanCrippenLogP import WildmanCrippenLogP as LogP, WildmanCrippenMR as MR
from rdkit import Chem
from numpy.testing import assert_almost_equal


def test_WildmanCrippen1():
    mol = Chem.MolFromSmiles('Oc1ccccc1OC')
    yield assert_almost_equal, LogP()(mol), 1.4, 2
    yield assert_almost_equal, MR()(mol), 34.66, 2


def test_WildmanCrippen2():
    mol = Chem.MolFromSmiles('c1ccccc1c2ccccn2')
    yield assert_almost_equal, LogP()(mol), 2.75, 2
