import math
import pickle

import six
from rdkit import Chem
from numpy.testing import assert_almost_equal

from mordred import Calculator, descriptors
from mordred.error import MissingValueBase


def test_pickle_calculator():
    orig = Calculator(descriptors)
    d0 = orig.descriptors[0]
    d1 = orig.descriptors[1]
    orig.register(
        [
            d0 + d1,
            d0 - d1,
            d0 * d1,
            d0 // d1,
            d0 % d1,
            d0**d1,
            -d0,
            +d1,
            abs(d0),
            math.trunc(d0),
        ]
    )

    if six.PY3:
        orig.register([math.ceil(d0), math.floor(d1)])

    pickled = pickle.loads(pickle.dumps(orig))

    mol = Chem.MolFromSmiles("c1ccccc1C(O)O")

    for a, b in zip(orig.descriptors, pickled.descriptors):
        assert a == b

    for a, b in zip(orig(mol), pickled(mol)):
        if isinstance(a, MissingValueBase):
            assert a.__class__ == b.__class__
        else:
            assert_almost_equal(a, b)
