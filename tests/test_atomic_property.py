from rdkit import Chem
from mold import _atomic_property

from numpy.testing import assert_almost_equal

_data = {
    '>C<': 'CC(C)(C)C',
    '>CH-': 'CC(C)C',
    '-CH2-': 'CCC',
    '=C<': 'C=C(C)C',
    '-CH3': 'CC',
    '=CH-': 'CC=C',
    '>N-': 'CN(C)C',
    '#C-': 'C#CC',
    '-NH-': 'CNC',
    '=CH2': 'C=C',
    '=N-': 'C=NC',
    '-O-': 'COC',
    '#CH': 'C#C',
    '-NH2': 'CN',
    '=NH': 'C=N',
    '#N': 'C#N',
    '-OH': 'CO',
    '=O': 'C=O',
    '-F': 'CF',
    '-SH': 'CS',
    '-S-': 'CSC',
    '=S': 'C=S',
    '-Cl': 'CCl',
    '-Br': 'CBr',
    '-I': 'CI',
}

explicitHs = {key: Chem.AddHs(Chem.MolFromSmiles(smi)).GetAtomWithIdx(1)
              for key, smi in _data.items()}

implicitHs = {key: Chem.MolFromSmiles(smi).GetAtomWithIdx(1)
              for key, smi in _data.items()}


def make_cases(results, fn, decimal=7):
    for key, result in results.items():
        yield assert_almost_equal, result, fn(explicitHs[key]), decimal, key
        yield assert_almost_equal, result, fn(implicitHs[key]), decimal, key


def test_sigma():
    results = {
        '>C<': 4,
        '>CH-': 3,
        '-CH2-': 2,
        '=C<': 3,
        '-CH3': 1,
        '=CH-': 2,
        '>N-': 3,
        '#C-': 2,
        '-NH-': 2,
        '=CH2': 1,
        '=N-': 2,
        '-O-': 2,
        '#CH': 1,
        '-NH2': 1,
        '=NH': 1,
        '#N': 1,
        '-OH': 1,
        '=O': 1,
        '-F': 1,
        '-SH': 1,
        '-S-': 2,
        '=S': 1,
        '-Cl': 1,
        '-Br': 1,
        '-I': 1,
    }
    for test in make_cases(results, _atomic_property.GetNumSigmaElectrons):
        yield test


def test_valence_sigma():
    results = {
        '>C<': 4,
        '>CH-': 3,
        '-CH2-': 2,
        '=C<': 4,
        '-CH3': 1,
        '=CH-': 3,
        '>N-': 5,
        '#C-': 4,
        '-NH-': 4,
        '=CH2': 2,
        '=N-': 5,
        '-O-': 6,
        '#CH': 3,
        '-NH2': 3,
        '=NH': 4,
        '#N': 5,
        '-OH': 5,
        '=O': 6,
        '-F': 7,
        '-S-': 0.67,
        '-Cl': 0.78,
        '-Br': 0.26,
        '-I': 0.16,
    }
    for test in make_cases(results, _atomic_property.GetNumValenceElectrons, 2):
        yield test


def test_IntrinsicState():
    results = {
        '>C<': 1.25,
        '>CH-': 1.3333,
        '-CH2-': 1.5,
        '=C<': 1.6666,
        '-CH3': 2.0,
        '=CH-': 2.0,
        '>N-': 2.0,
        '#C-': 2.5,
        '-NH-': 2.5,
        '=CH2': 3.0,
        '=N-': 3.0,
        '-O-': 3.5,
        '#CH': 4.0,
        '-NH2': 4.0,
        '=NH': 5.0,
        '#N': 6.0,
        '-OH': 6.0,
        '=O': 7.0,
        '-F': 8.0,
    }

    for test in make_cases(results, _atomic_property.IntrinsicState, 3):
        yield test
