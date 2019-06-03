import math
import operator

import six
from rdkit import Chem

from nose.tools import eq_

from .. import ABCIndex

binary = [
    operator.add,
    operator.sub,
    operator.mul,
    operator.truediv,
    operator.floordiv,
    operator.mod,
    operator.pow,
]

unary = [operator.neg, operator.pos, operator.abs, math.trunc]


if six.PY3:
    unary.extend([math.ceil, math.floor])


def test_compose_descriptor():
    L = ABCIndex.ABCIndex()
    r = ABCIndex.ABCGGIndex()
    m = Chem.MolFromSmiles("c1ccccc1C")

    for op in binary:
        yield eq_, op(L, r)(m), op(L(m), r(m))

    for op in unary:
        yield eq_, op(L)(m), op(L(m))
