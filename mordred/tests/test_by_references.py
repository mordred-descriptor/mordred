import os
from glob import glob

import yaml
import numpy as np
from rdkit import Chem
from numpy.testing import assert_almost_equal

from mordred import Calculator, Polarizability, descriptors
from mordred.error import MissingValueBase

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


data_dir = os.path.join(os.path.dirname(__file__), "references")


def test_by_references():
    calc = Calculator(
        d
        for d in descriptors.all
        if d.__class__ not in [Polarizability.APol, Polarizability.BPol]
    )

    calc.register([Polarizability.APol(True), Polarizability.BPol(True)])

    actuals = {}
    for mol in Chem.SDMolSupplier(
        os.path.join(data_dir, "structures.sdf"), removeHs=False
    ):
        actuals[mol.GetProp("_Name")] = {
            str(d): v for d, v in zip(calc.descriptors, calc(mol))
        }

    for path in glob(os.path.join(data_dir, "*.yaml")) + glob(
        os.path.join(data_dir, "**/*.yaml")
    ):
        for test in yaml.load(open(path), Loader=Loader):
            dnames = test["names"]
            if not isinstance(dnames, list):
                dnames = [dnames]

            desireds = (
                (mname, zip(dnames, values if isinstance(values, list) else [values]))
                for mname, values in test["results"].items()
            )

            digit = test.get("digit")
            if digit is None:
                digit = 0  # exact match

            def assert_f(a, d, m):
                if np.isnan(d):
                    assert isinstance(a, MissingValueBase)
                    return

                assert_almost_equal(a, d, digit, m)

            for mname, descs in desireds:
                for dname, desired in descs:
                    if not desired == "skip":
                        assert_f(
                            actuals[mname][dname],
                            desired,
                            "{} of {}".format(dname, mname),
                        )
