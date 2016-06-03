from rdkit import Chem
from mordred import Descriptor, Calculator
import os

from nose.tools import eq_


class DummyLogger(object):
    def warning(self, *args):
        eq_(args[1], os.path.basename(__file__))
        eq_(args[3], 'test exception')


class ExceptionDescriptor(Descriptor):
    def __reduce_ex__(self, version):
        return self.__class__, ()

    def calculate(self, mol):
        raise ValueError('test exception')


mol = Chem.MolFromSmiles('c1ccccc1')


def test_catch_exception():
    calc = Calculator(ExceptionDescriptor())
    calc.logger = DummyLogger()
    calc(mol)
