import argparse
import csv
import sys

from rdkit import Chem

from . import Calculator, all_descriptors


def smiles_parser(f):
    for line in f:
        smi, name = line.strip().split()
        mol = Chem.MolFromSmiles(smi)
        mol.SetProp('_Name', name)
        yield mol


def main(descs):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', metavar='PATH',
        nargs='?', type=argparse.FileType('r'), default=sys.stdin,
        help='input smiles file(default: stdin)',
    )

    parser.add_argument(
        '-o', '--output', metavar='PATH',
        nargs='?', type=argparse.FileType('w'), default=sys.stdout,
        help='output csv file(default: stdout)'
    )

    parser.add_argument(
        '-p', '--processes', metavar='N',
        default=None, type=int,
        help='number of processes to use',
    )

    args = parser.parse_args()

    mols = list(smiles_parser(args.input))

    calc = Calculator(descs)

    with args.output as output:
        writer = csv.writer(output)

        writer.writerow(['name'] + list(map(str, calc.descriptors)))

        for mol, val in zip(mols, calc.map(mols, args.processes)):
            writer.writerow([mol.GetProp('_Name')] + list(map(lambda r: str(r[1]), val)))

    return 0


def submodule(descs):
    sys.exit(main(descs))


if __name__ == '__main__':
    sys.exit(main(all_descriptors()))
