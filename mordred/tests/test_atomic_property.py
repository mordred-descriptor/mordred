from rdkit import Chem
from numpy.testing import assert_almost_equal

from mordred import _atomic_property

_data = {
    ">C<": "CC(C)(C)C",
    ">CH-": "CC(C)C",
    "-CH2-": "CCC",  # noqa: S001
    "=C<": "C=C(C)C",
    "-CH3": "CC",  # noqa: S001
    "=CH-": "CC=C",
    ">N-": "CN(C)C",
    "#C-": "C#CC",  # noqa: S001
    "-NH-": "CNC",
    "=CH2": "C=C",
    "=N-": "C=NC",
    "-O-": "COC",  # noqa: S001
    "#CH": "C#C",  # noqa: S001
    "-NH2": "CN",
    "=NH": "C=N",
    "#N": "C#N",  # noqa: S001
    "-OH": "CO",
    "=O": "C=O",
    "-F": "CF",  # noqa: S001
    "-SH": "CS",
    "-S-": "CSC",  # noqa: S001
    "=S": "C=S",
    "-Cl": "CCl",  # noqa: S001
    "-Br": "CBr",  # noqa: S001
    "-I": "CI",
}

explicitHs = {
    key: Chem.AddHs(Chem.MolFromSmiles(smi)).GetAtomWithIdx(1)
    for key, smi in _data.items()
}

implicitHs = {
    key: Chem.MolFromSmiles(smi).GetAtomWithIdx(1) for key, smi in _data.items()
}


def make_cases(results, fn, decimal=7):
    for key, result in results.items():
        assert assert_almost_equal(result, fn(explicitHs[key]), decimal, key) is None
        assert assert_almost_equal(result, fn(implicitHs[key]), decimal, key) is None


def test_sigma():
    results = {
        ">C<": 4,
        ">CH-": 3,
        "-CH2-": 2,  # noqa: S001
        "=C<": 3,
        "-CH3": 1,  # noqa: S001
        "=CH-": 2,
        ">N-": 3,
        "#C-": 2,  # noqa: S001
        "-NH-": 2,
        "=CH2": 1,
        "=N-": 2,
        "-O-": 2,  # noqa: S001
        "#CH": 1,  # noqa: S001
        "-NH2": 1,
        "=NH": 1,
        "#N": 1,  # noqa: S001
        "-OH": 1,
        "=O": 1,
        "-F": 1,  # noqa: S001
        "-SH": 1,
        "-S-": 2,  # noqa: S001
        "=S": 1,
        "-Cl": 1,  # noqa: S001
        "-Br": 1,  # noqa: S001
        "-I": 1,
    }
    make_cases(results, _atomic_property.get_sigma_electrons)


def test_valence_sigma():
    results = {
        ">C<": 4,
        ">CH-": 3,
        "-CH2-": 2,  # noqa: S001
        "=C<": 4,
        "-CH3": 1,  # noqa: S001
        "=CH-": 3,
        ">N-": 5,
        "#C-": 4,  # noqa: S001
        "-NH-": 4,
        "=CH2": 2,
        "=N-": 5,
        "-O-": 6,  # noqa: S001
        "#CH": 3,  # noqa: S001
        "-NH2": 3,
        "=NH": 4,
        "#N": 5,  # noqa: S001
        "-OH": 5,
        "=O": 6,
        "-F": 7,  # noqa: S001
        "-S-": 0.67,
        "-Cl": 0.78,  # noqa: S001
        "-Br": 0.26,  # noqa: S001
        "-I": 0.16,
    }
    make_cases(results, _atomic_property.get_valence_electrons, 2)


def test_IntrinsicState():
    results = {
        ">C<": 1.25,
        ">CH-": 1.3333,
        "-CH2-": 1.5,  # noqa: S001
        "=C<": 1.6666,
        "-CH3": 2.0,  # noqa: S001
        "=CH-": 2.0,
        ">N-": 2.0,
        "#C-": 2.5,  # noqa: S001
        "-NH-": 2.5,
        "=CH2": 3.0,
        "=N-": 3.0,
        "-O-": 3.5,  # noqa: S001
        "#CH": 4.0,  # noqa: S001
        "-NH2": 4.0,
        "=NH": 5.0,
        "#N": 6.0,  # noqa: S001
        "-OH": 6.0,
        "=O": 7.0,
        "-F": 8.0,  # noqa: S001
    }

    make_cases(results, _atomic_property.get_intrinsic_state, 3)
