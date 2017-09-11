from rdkit import Chem

from mordred import Calculator, descriptors

mol = Chem.MolFromSmiles("c1ccccc1")

# Start Code 7
calc = Calculator(descriptors)
result = calc(mol)
result_dict = result.drop_missing().asdict()
print(len(result_dict))
print(result_dict["SLogP"])
