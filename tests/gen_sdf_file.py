import os

from rdkit.Chem import AllChem as Chem

base = os.path.join(
    os.path.dirname(__file__),
    'references', 'structures'
)

writer = Chem.SDWriter(base + '.sdf')
for line in open(base + '.smi'):
    smi, name = line.strip().split()
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    mol.SetProp('_Name', name)
    
    Chem.EmbedMolecule(mol, randomSeed=23216)
    while Chem.MMFFOptimizeMolecule(mol) != 0:
        pass

    writer.write(mol)
    writer.flush()
