import json

from mordred import Calculator, descriptors
from nose.tools import eq_


def test_json():
    calc = Calculator(descriptors)
    j = json.dumps(calc.to_json())
    calc2 = Calculator.from_json(json.loads(j))
    eq_(calc.descriptors, calc2.descriptors)
