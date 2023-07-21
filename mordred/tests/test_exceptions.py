from rdkit import Chem

from mordred import Calculator, Descriptor
from mordred.error import Error


class RaiseDescriptor(Descriptor):
    def parameters(self):
        return ()

    def __init__(self, e, critical):
        self.e = e
        self.critical = critical

    def calculate(self):
        raise self.e


mol = Chem.MolFromSmiles("c1ccccc1")


def test_catch_non_critical_error():
    calc = Calculator(RaiseDescriptor(ValueError("test exception"), False))

    result = calc(mol)[0]
    assert isinstance(result, Error)
    assert result.error.args == ("test exception",)
