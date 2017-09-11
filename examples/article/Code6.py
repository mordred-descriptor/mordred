from rdkit import Chem

from mordred import Calculator
from mordred.GeometricalIndex import Radius3D

mol2D = Chem.MolFromSmiles("c1ccccc1")

# Start Code 6
calc = Calculator(Radius3D)
err, = calc(mol2D)
print(repr(err))
print(err.error)
