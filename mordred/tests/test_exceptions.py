from rdkit import Chem
from mordred import Descriptor, Calculator
from mordred.exception import MordredException
import os

from nose.tools import eq_, raises


class DummyLogger(object):
    def warning(self, *args):
        eq_(args[1], os.path.basename(__file__))
        eq_(args[3], 'test exception')


class CriticalError(MordredException):
    critical = True


class RaiseDescriptor(Descriptor):
    def as_key(self):
        return self.__class__, ()

    def __init__(self, e):
        self.e = e

    def calculate(self, mol):
        raise self.e


mol = Chem.MolFromSmiles('c1ccccc1')


@raises(ValueError)
def test_through_unknown_exception():
    calc = Calculator(RaiseDescriptor(ValueError()))
    calc(mol)


def test_catch_non_critical_error():
    calc = Calculator(RaiseDescriptor(MordredException('test exception')))
    calc.logger = DummyLogger()
    calc(mol)


@raises(CriticalError)
def test_through_critical_error():
    calc = Calculator(RaiseDescriptor(CriticalError()))
    calc(mol)
