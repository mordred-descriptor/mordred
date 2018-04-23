# Based on the code from https://github.com/rdkit/rdkit/blob/master/Contrib/PBF/pbf.py
from numpy.testing import assert_almost_equal
from rdkit import Chem

from mordred.SP3A2C import SP3A2C


def test_SP3A2C():
    m = Chem.MolFromSmiles('O=C(CC#N)N1CCC[C@H](C1)Nc1ccnc(n1)-c1cnc2ccc(cn12)C#N')
    yield assert_almost_equal, SP3A2C()(m), 0.833, 2
