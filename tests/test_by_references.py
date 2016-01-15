import os
from glob import glob
from nose.tools import eq_
from numpy.testing import assert_almost_equal
from rdkit import Chem

from mordred import Calculator, all_descriptors
from mordred import Polarizability

import yaml
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


data_dir = os.path.join(
    os.path.dirname(__file__),
    'data'
)


def test_by_references():
    calc = Calculator(all_descriptors())

    calc.descriptors = list(
        filter(lambda x: x.__class__ not in [Polarizability.APol,
                                             Polarizability.BPol], calc.descriptors))
    calc.register(Polarizability.APol(True),
                  Polarizability.BPol(True),
                  )

    actuals = dict()
    for mol in Chem.SmilesMolSupplier(os.path.join(data_dir, 'test.smi'), titleLine=False):
        actuals[mol.GetProp('_Name')] = dict(calc(mol))

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
                assert_f = lambda a, d, m: assert_almost_equal(a, d, digit, m)

            for mname, descs in desireds:
                for dname, desired in descs:
                    if not desired == 'skip':
                        yield assert_f,\
                            actuals[mname][dname],\
                            desired,\
                            '{} of {}'.format(dname, mname)
