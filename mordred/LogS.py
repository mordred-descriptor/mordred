from ._base import Descriptor
from .Weight import Weight
import math
import rdkit.Chem as Chem

__all__ = (
    "LogS",
)

smarts_logs = {
    '[NH0;X3;v3]': 0.71535,
    '[NH2;X3;v3]': 0.41056,
    '[nH0;X3]': 0.82535,
    '[OH0;X2;v2]': 0.31464,
    '[OH0;X1;v2]': 0.14787,
    '[OH1;X2;v2]': 0.62998,
    '[CH2;!R]': -0.35634,
    '[CH3;!R]': -0.33888,
    '[CH0;R]': -0.21912,
    '[CH2;R]': -0.23057,
    '[ch0]': -0.37570,
    '[ch1]': -0.22435,
    'F': -0.21728,
    'Cl': -0.49721,
    'Br': -0.57982,
    'I': -0.51547
}


class LogS(Descriptor):
    r"""LogS descriptor.
    """

    __slots__ = ("_type",)

    @classmethod
    def preset(cls):
        yield cls()

    def dependencies(self):
        return {
            "MW": Weight(),
        }

    def description(self):
        return "LogS"

    def __str__(self):
        return self._type

    def __init__(self, type="LogS"):
        self._type = type

    def parameters(self):
        return ()

    def calculate(self, MW):
        logS = 0.89823 - 0.10369 * math.sqrt(MW)

        for key in smarts_logs:
            logS += len(self.mol.GetSubstructMatches(Chem.MolFromSmarts(key))) * smarts_logs[key]

        return logS

    rtype = int
