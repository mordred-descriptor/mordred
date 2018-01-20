from rdkit import Chem

from mordred import Calculator
from mordred.GeometricalIndex import Radius3D

# Start Code 6
calc = Calculator(Radius3D)
result = calc(Chem.MolFromSmiles("c1ccccc1"))
err = result[0]
print(repr(err))
print(err.error)
