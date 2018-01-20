from rdkit import Chem

from mordred.SLogP import SLogP
from mordred.Lipinski import Lipinski

# Start Code 8
slogp = SLogP()
lipinski = Lipinski()
product_term = slogp * lipinski
print(str(product_term))
print(product_term(Chem.MolFromSmiles("c1ccccc1")))
