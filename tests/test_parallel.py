from mordred import Calculator, all_descriptors
from rdkit import Chem
from nose.tools import eq_
from numpy.testing import assert_almost_equal


def test_parallel():
    calc = Calculator(all_descriptors())
    mol = Chem.MolFromSmiles('C1=CC(=C(C=C1C2=C(C=C3C(=CC(=CC3=[O+]2)O)O)O)O)O')

    serial = calc(mol)
    parallel = calc.parallel([mol])[0]

    eq_([kv[0] for kv in serial], [kv[0] for kv in parallel])

    for (_, s), (_, p) in zip(serial, parallel):
        assert_almost_equal(s, p)
