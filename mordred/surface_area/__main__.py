from __future__ import print_function

import argparse

from rdkit import Chem

from . import SurfaceArea
from .._util import PathType, module_prog


def main():
    parser = argparse.ArgumentParser(
        prog=module_prog(__package__), formatter_class=argparse.MetavarTypeHelpFormatter
    )
    parser.add_argument("sdf", type=PathType, help="input sd file", metavar="SDF")
    parser.add_argument(
        "-s",
        "--solvent-radius",
        type=float,
        default=1.4,
        help="solvent radius (default: %(default)s)",
    )
    parser.add_argument(
        "-l",
        "--mesh-level",
        type=int,
        default=5,
        help="mesh level (default: %(default)s)",
    )
    result = parser.parse_args()
    main_process(
        sdf=result.sdf,
        solvent_radius=result.solvent_radius,
        mesh_level=result.mesh_level,
    )


def main_process(sdf, solvent_radius, mesh_level):
    for i, mol in enumerate(Chem.SDMolSupplier(sdf, removeHs=False)):
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else str(i)
        sa = SurfaceArea.from_mol(mol, solvent_radius=solvent_radius, level=mesh_level)

        print(name)  # noqa: T003

        for a, s in zip(mol.GetAtoms(), sa.surface_area()):
            print(  # noqa: T001
                "{:4d} {:2s} {:8.3f}".format(a.GetIdx(), a.GetSymbol(), s)
            )  # noqa: T003


if __name__ == "__main__":
    main()
