from mordred import Calculator, all_descriptors
from rdkit import Chem
from nose.tools import eq_
from numpy.testing import assert_almost_equal


def test_parallel():
    calc = Calculator(all_descriptors())
    mols = [
        Chem.MolFromSmiles('c1ccccc1'),
        Chem.MolFromSmiles('C1=CC(=C(C=C1C2=C(C=C3C(=CC(=CC3=[O+]2)O)O)O)O)O'),
        Chem.MolFromSmiles('CCCCCC'),
    ]

    for (_, serial), (_, parallel) in zip(calc.map(mols, processes=1), calc.map(mols)):
        for (d, s), (f, p) in zip(serial, parallel):
            yield eq_, d, f
            yield assert_almost_equal, s, p
