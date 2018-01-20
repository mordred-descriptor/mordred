from rdkit import Chem

from mordred.RingCount import RingCount

# Start Code 2
rings = RingCount()
rings6 = RingCount(order=6)
print(rings(Chem.MolFromSmiles("c1ccccc1")), rings(Chem.MolFromSmiles("C1CCCCCC1")))
print(rings6(Chem.MolFromSmiles("c1ccccc1")), rings6(Chem.MolFromSmiles("C1CCCCCC1")))
