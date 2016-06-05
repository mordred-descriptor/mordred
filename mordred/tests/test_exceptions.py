from rdkit import Chem
from mordred import Descriptor, Calculator, Error

from nose.tools import raises, eq_


class RaiseDescriptor(Descriptor):
    def as_key(self):
        return self.__class__, ()

    def __init__(self, e, critical):
        self.e = e
        self.critical = critical

    def calculate(self):
        self.fail(self.e, critical=self.critical)


mol = Chem.MolFromSmiles('c1ccccc1')


@raises(Error)
def test_through_critical_error():
    calc = Calculator(RaiseDescriptor(ValueError(), True))
    calc(mol)


def test_catch_non_critical_error():
    calc = Calculator(RaiseDescriptor(ValueError('test exception'), False))

    result = calc(mol)[0]
    assert isinstance(result, Error)
    eq_(result.error.args, ('test exception',))
