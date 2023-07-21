import os

from rdkit.Chem import AllChem as Chem
from numpy.testing import assert_almost_equal

from mordred import Calculator, descriptors
from mordred.error import MissingValueBase

data_file = os.path.join(os.path.dirname(__file__), "references", "structures.sdf")


def test_parallel():
    calc = Calculator(descriptors)
    mols = [m for m in Chem.SDMolSupplier(data_file, removeHs=False)]

    for serial, parallel in zip(
        calc.map(mols, nproc=1, quiet=True), calc.map(mols, quiet=True)
    ):
        for d, s, p in zip(calc.descriptors, serial, parallel):
            if isinstance(s, MissingValueBase):
                assert s.error.__class__ == p.error.__class__
            else:
                msg = "{} (serial: {}, parallel: {})".format(str(d), s, p)
                assert assert_almost_equal(s, p, 7, msg) is None
