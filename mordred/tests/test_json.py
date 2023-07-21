import json

from mordred import Calculator, descriptors


def test_json():
    calc = Calculator(descriptors)
    j = json.dumps(calc.to_json())
    calc2 = Calculator.from_json(json.loads(j))
    assert calc.descriptors == calc2.descriptors
