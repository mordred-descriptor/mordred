import os
from glob import glob
from nose.tools import eq_
from numpy.testing import assert_almost_equal
from rdkit import Chem

import numpy as np
from mordred import Calculator, all_descriptors, Polarizability
from mordred.error import MissingValue

import yaml
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


data_dir = os.path.join(
    os.path.dirname(__file__),
    'references'
)


def test_by_references():
    calc = Calculator(
        d
        for d in all_descriptors()
        if d.__class__ not in [Polarizability.APol, Polarizability.BPol]
    )

    calc.register([
        Polarizability.APol(True),
        Polarizability.BPol(True),
    ])

    actuals = dict()
    for mol in Chem.SDMolSupplier(os.path.join(data_dir, 'structures.sdf'), removeHs=False):
        actuals[mol.GetProp('_Name')] = {str(d): v for d, v in zip(calc.descriptors, calc(mol))}

    for path in glob(os.path.join(data_dir, '*.yaml')) + glob(os.path.join(data_dir, '**/*.yaml')):
        for test in yaml.load(open(path), Loader=Loader):
            dnames = test['names']
            if not isinstance(dnames, list):
                dnames = [dnames]

            desireds = (
                (mname, zip(dnames, values if isinstance(values, list) else [values]))
                for mname, values in test['results'].items()
            )

            digit = test.get('digit')
            if digit is None:
                assert_f = eq_
            else:
                def assert_f(a, d, m):
                    if np.isnan(d):
                        assert isinstance(a, MissingValue)
                        return

                    assert_almost_equal(a, d, digit, m)

            for mname, descs in desireds:
                for dname, desired in descs:
                    if not desired == 'skip':
                        yield (
                            assert_f,
                            actuals[mname][dname],
                            desired,
                            '{} of {}'.format(dname, mname)
                        )
