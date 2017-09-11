from rdkit import Chem

from mordred.RingCount import RingCount

# Start Code 3
rings50 = RingCount(order=50)
print(rings50(Chem.MolFromSmiles("c1ccccc1")), rings50(Chem.MolFromSmiles("C1CCCCCC1")))
