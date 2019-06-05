from rdkit import Chem

from mordred import Chi, ABCIndex

benzene = Chem.MolFromSmiles('c1ccccc1')

# create descriptor instance
abci = ABCIndex.ABCIndex()

# calculate descriptor value
result = abci(benzene)

print(str(abci), result)

# create descriptor instance with parameter
chi_pc4 = Chi.Chi(type='path_cluster', order=4)

# calculate
result = chi_pc4(benzene)

print(str(chi_pc4), result)
