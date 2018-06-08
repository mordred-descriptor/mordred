from rdkit import Chem
from . import get_molecular_framework
from itertools import chain


def main():
    import sys
    mol = Chem.MolFromSmiles(sys.argv[1])
    linkers, rings = get_molecular_framework(mol)
    indices = set(chain(chain(*linkers), chain(*rings)))

    rwmol = Chem.EditableMol(mol)
    for i in range(mol.GetNumAtoms() - 1, -1, -1):
        if i not in indices:
            rwmol.RemoveAtom(i)

    print(Chem.MolToSmiles(rwmol.GetMol()))


if __name__ == "__main__":
    main()
