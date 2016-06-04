from rdkit import Chem
from mordred import Descriptor, Calculator
from mordred.exception import MordredException

from nose.tools import raises, eq_


class DummyLogger(object):
    def __init__(self):
        self.results = []

    def warning(self, *args):
        self.results.append(args)


class CriticalError(MordredException):
    critical = True


class RaiseDescriptor(Descriptor):
    def as_key(self):
        return self.__class__, ()

    def __init__(self, e):
        self.e = e

    def calculate(self):
        raise self.e


mol = Chem.MolFromSmiles('c1ccccc1')


@raises(ValueError)
def test_through_unknown_exception():
    calc = Calculator(RaiseDescriptor(ValueError()))
    calc(mol)


def test_catch_non_critical_error():
    calc = Calculator(RaiseDescriptor(MordredException('test exception')))
    calc.logger = DummyLogger()

    result = calc(mol)[0]
    assert isinstance(result, MordredException)
    eq_(result.args, ('test exception',))


@raises(CriticalError)
def test_through_critical_error():
    calc = Calculator(RaiseDescriptor(CriticalError()))
    calc(mol)
