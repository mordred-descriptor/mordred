import math
import operator

import six
from rdkit import Chem


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
        assert op(L, r)(m) == op(L(m), r(m))

    for op in unary:
        assert op(L)(m) == op(L(m))
