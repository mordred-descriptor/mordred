from rdkit import Chem
from numpy.testing import assert_almost_equal

from mordred.SLogP import SMR, SLogP


def test_WildmanCrippen1():
    mol = Chem.MolFromSmiles("Oc1ccccc1OC")
    assert assert_almost_equal(SLogP()(mol), 1.4, 2) is None
    assert assert_almost_equal(SMR()(mol), 34.66, 2) is None


def test_WildmanCrippen2():
    mol = Chem.MolFromSmiles("c1ccccc1c2ccccn2")
    assert assert_almost_equal(SLogP()(mol), 2.75, 2) is None
