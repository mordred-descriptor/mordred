from rdkit import Chem

from mordred import RingCount, Calculator

desc = RingCount.RingCount()
mol = Chem.MolFromSmiles("c1ccccc1")
mols = [
    Chem.MolFromSmiles(smi)
    for smi in [
        "CCCCCC",
        "C1CCCC1",
        "C1CCCCC1",
    ]
]

# Start Code 5
calc = Calculator()
calc.register(desc)
print(calc(mol))
print(list(calc.map(mols)))
