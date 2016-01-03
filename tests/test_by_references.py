import os
from mold import Calculator
import mold.all as md
from functools import wraps

from rdkit import Chem
from math import isnan
import numpy.testing as nt

import yaml
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


ref_file = os.path.join(
    os.path.dirname(__file__),
    'data', 'test.yaml'
)


def assert_equal(actual, desired, name, desc):
    assert actual == desired,\
        '{!r} != {!r}, {} of {}'.format(actual, desired, desc, name)


def assert_almost_equal(decimal):
    def f(actual, desired, name, desc):
        nt.assert_almost_equal(
            actual, desired, decimal, '{} of {}'.format(desc, name))
    return f


def fill_na(v, f):
    @wraps(f)
    def g(actual, desired, name, desc):
        actual = v if isnan(actual) else actual
        desired = v if isnan(desired) else desired

        f(actual, desired, name, desc)

    return g


def all_(a):
    return True


def exclude_(*ns):
    def f(name):
        return name not in ns

    return f


def only_(*ns):
    def f(name):
        return name in ns

    return f

descriptors = [
    (assert_equal, md.Aromatic, all_),
    (assert_equal, md.AtomCount, all_),
    (assert_almost_equal(7),
     [[md.Autocorrelation.ATS(d, a),
       md.Autocorrelation.ATSC(d, a),
       ] for a in 'mvepis' for d in range(9)], exclude_('Cyanidin')
     ),
    (fill_na(0, assert_almost_equal(7)),
     [[md.Autocorrelation.AATS(d, a),
       md.Autocorrelation.AATSC(d, a),
       ] for a in 'mvepis' for d in range(9)], exclude_('Cyanidin')
     ),
    (fill_na(0, assert_almost_equal(7)),
     [[md.Autocorrelation.MATS(d, a),
       md.Autocorrelation.GATS(d, a),
       ] for a in 'mvepis' for d in range(1, 9)], exclude_('Cyanidin')
     ),

    # (assert_equal, md.BondCount, []),

    (assert_equal, md.CarbonTypes, all_),

    (assert_almost_equal(0),
     [md.Matrix.BCUT('m', 1, False),
      md.Matrix.BCUT('m', 1, True),
      ], exclude_('Cyanidin')),

    (assert_almost_equal(7),
     [md.Matrix.BaryszMatrix('Z', a)
      for a in ['SpAbs', 'SpMax', 'SpDiam',
                'SpAD', 'SpMAD',
                'EE', 'SM1',
                'VE1', 'VE2',
                'VR1', 'VR2', 'VR3']
      ], only_('Hexane')),

    (assert_almost_equal(7),
     [md.Chi.Chi('path', l, a) for a in ['delta', 'delta_v'] for l in range(7)] +
     [md.Chi.Chi('chain', l, a) for a in ['delta', 'delta_v'] for l in range(3, 7)] +
     [md.Chi.Chi('cluster', l, a) for a in ['delta', 'delta_v'] for l in range(3, 7)] +
     [md.Chi.Chi('path_cluster', l, a) for a in ['delta', 'delta_v'] for l in range(4, 6)],
     exclude_('Cyanidin')),

    (assert_almost_equal(7), [md.Polarizability.APol(True),
                              md.Polarizability.BPol(True),
                              ], exclude_('Cyanidin')),

    (assert_equal, md.SmartsCount, all_),
]


def test_by_references():
    for data in yaml.load(open(ref_file), Loader=Loader):
        name = data['name']
        smi = data['smiles']
        descs = data['descriptors']

        for check, desc, is_check in descriptors:
            if not is_check(name):
                continue

            calc = Calculator(desc)
            for desc_name, actual in calc(Chem.MolFromSmiles(smi)):
                desired = descs[desc_name]

                yield check, actual, desired, name, desc_name
