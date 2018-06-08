from logging import StreamHandler, getLogger

from rdkit import Chem

from . import get_molecular_framework

logger = getLogger(__file__)


def main():
    import sys
    mol = Chem.MolFromSmiles(sys.argv[1])
    linkers, rings = get_molecular_framework(mol)
    indices = {i for linker in linkers for ab in linker for i in ab}
    indices.update(i for ring in rings for i in ring)

    rwmol = Chem.EditableMol(mol)
    for i in range(mol.GetNumAtoms() - 1, -1, -1):
        if i not in indices:
            rwmol.RemoveAtom(i)

    logger.info(Chem.MolToSmiles(rwmol.GetMol()))


if __name__ == "__main__":
    logger.addHandler(StreamHandler())
    logger.setLevel(20)
    main()
