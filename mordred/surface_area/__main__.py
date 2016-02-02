import argparse

from rdkit import Chem

from . import SurfaceArea


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'FILE', type=str,
        help='input sdf/mol file',
    )

    parser.add_argument(
        '-s', '--solvent-radius',
        type=float, default=1.4,
        help='solvent radius (default: %(default)s)',
    )

    parser.add_argument(
        '-l', '--mesh-level',
        type=int, default=5,
        help='mesh level (default: %(default)s)',
    )

    args = parser.parse_args()

    for i, mol in enumerate(Chem.SDMolSupplier(args.FILE, removeHs=False)):
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else str(i)
        sa = SurfaceArea.from_mol(mol, solvent_radius=args.solvent_radius, level=args.mesh_level)

        print(name)

        for a, s in zip(mol.GetAtoms(), sa.surface_area()):
            print('{:4d} {:2s} {:8.3f}'.format(
                a.GetIdx(),
                a.GetSymbol(),
                s
            ))


if __name__ == '__main__':
    main()
