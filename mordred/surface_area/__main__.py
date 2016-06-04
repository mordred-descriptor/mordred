import sys
import os
import click

from rdkit import Chem

from . import SurfaceArea


@click.command(
    context_settings={
        'help_option_names': ['-h', '--help']
    }
)
@click.argument('SDF', type=click.Path(exists=True))
@click.option(
    '-s', '--solvent-radius',
    type=click.FLOAT, default=1.4,
    help='solvent radius'
)
@click.option(
    '-l', '--mesh-level',
    type=click.INT, default=5,
    help='mesh level'
)
def main(sdf, solvent_radius, mesh_level):
    for i, mol in enumerate(Chem.SDMolSupplier(sdf, removeHs=False)):
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else str(i)
        sa = SurfaceArea.from_mol(mol, solvent_radius=solvent_radius, level=mesh_level)

        print(name)

        for a, s in zip(mol.GetAtoms(), sa.surface_area()):
            print('{:4d} {:2s} {:8.3f}'.format(
                a.GetIdx(),
                a.GetSymbol(),
                s
            ))


if __name__ == '__main__':
    main(prog_name='{} -m {}'.format(os.path.basename(sys.executable), __package__))
