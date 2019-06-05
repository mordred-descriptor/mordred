from rdkit import Chem

from mordred import Chi, ABCIndex

benzene = Chem.MolFromSmiles('c1ccccc1')

# create descriptor instance
abci = ABCIndex.ABCIndex()
chi_p2 = Chi.Chi(type='path', order=2)

# create product term using descriptor arithmetic
abci_x_chi_p2 = abci * chi_p2

# calculate descriptor value
result = abci_x_chi_p2(benzene)

print(abci_x_chi_p2, result)
