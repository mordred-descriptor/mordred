import os

from rdkit.Chem import AllChem as Chem

BASE = os.path.join(
    os.path.dirname(__file__),
    '..', 'mordred', 'tests', 'references'
)


def main():
    smi = os.path.join(BASE, 'structures.smi')
    sdf = os.path.join(BASE, 'structures.sdf')

    writer = Chem.SDWriter(sdf)
    for line in open(smi):
        smi, name = line.strip().split()
        print(name)
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        mol.SetProp('_Name', name)

        Chem.EmbedMolecule(mol, randomSeed=23216)
        while Chem.MMFFOptimizeMolecule(mol) == 1:
            pass

        writer.write(mol)
        writer.flush()


if __name__ == '__main__':
    main()
