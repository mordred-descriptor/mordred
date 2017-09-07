from multiprocessing import freeze_support

from rdkit import Chem

from mordred import Calculator, descriptors

freeze_support()

mols = [
    Chem.MolFromSmiles("c1ccccc1"),
    Chem.MolFromSmiles("c1ccccc1Cl"),
    Chem.MolFromSmiles("c1ccccc1C"),
]

# Create Calculator
calc = Calculator(descriptors)

# map method calculate multiple molecules (return generator)
print(list(calc.map(mols)))

# pandas method calculate multiple molecules (return pandas DataFrame)
print(calc.pandas(mols))
