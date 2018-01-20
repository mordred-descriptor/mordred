from rdkit import Chem

from mordred import Calculator, descriptors

# Start Code 7
calc = Calculator(descriptors)
result = calc(Chem.MolFromSmiles("c1ccccc1"))
result_dict = result.drop_missing().asdict()
print(len(result_dict))
print(result_dict["SLogP"])
