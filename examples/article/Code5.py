from multiprocessing import freeze_support

from rdkit import Chem

from mordred import Calculator
from mordred.RingCount import RingCount

# for parallel calculation on windows only
if __name__ == "__main__":
    freeze_support()

    # Start Code 5
    calc = Calculator()
    calc.register(RingCount())
    print(calc(Chem.MolFromSmiles("c1ccccc1")))
    print(list(calc.map([Chem.MolFromSmiles("c1ccccc1"), Chem.MolFromSmiles("CCCCCC")])))
