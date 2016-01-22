import argparse
import csv
import os
import sys

from rdkit import Chem

from ._base import Calculator, all_descriptors, get_descriptors_from_module


def smiles_parser(f):
    for line in f:
        line = line.strip().split()
        if len(line) == 1:
            smi = line[0]
            name = smi
        else:
            smi = line[0]
            name = ' '.join(line[1:])

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            sys.stderr.write('read failure: {}\n'.format(name))
            continue

        mol.SetProp('_Name', name)
        yield mol


def sdf_parser(p):
    for mol in Chem.SDMolSupplier(p, removeHs=False):
        yield mol


def file_parser(ifile, fmt):
    if fmt == 'smi':
        it = smiles_parser(ifile)
    elif fmt == 'sdf':
        ifile.close()
        it = sdf_parser(ifile.name)
    elif fmt == 'auto':
        ext = os.path.splitext(ifile.name)[1][1:]
        it = file_parser(ifile, ext)
    else:
        raise ValueError('unknown format: {!r}'.format(fmt))

    for mol in it:
        yield mol


def main(descs, prog=None):
    parser_options = dict()
    if prog is not None:
        parser_options['prog'] = prog

    parser = argparse.ArgumentParser(**parser_options)

    parser.add_argument(
        '-i', '--input', metavar='PATH',
        nargs='?', type=argparse.FileType('r'), default=sys.stdin,
        help='input file(default: stdin)',
    )

    parser.add_argument(
        '-f', '--from', metavar='TYPE',
        default='auto', choices=['auto', 'smi', 'sdf'],
        help='input filetype(one of %(choices)s, default: %(default)s)',
    )

    parser.add_argument(
        '-o', '--output', metavar='PATH',
        nargs='?', type=argparse.FileType('w'), default=sys.stdout,
        help='output csv file(default: stdout)'
    )

    parser.add_argument(
        '-p', '--processes', metavar='N',
        default=None, type=int,
        help='number of processes to use(default: number of threads)',
    )

    args = parser.parse_args()

    mols = file_parser(args.input, getattr(args, 'from'))

    calc = Calculator(descs)

    with args.output as output:
        writer = csv.writer(output)

        writer.writerow(['name'] + list(map(str, calc.descriptors)))

        for mol, val in calc.map(mols, args.processes):
            writer.writerow([mol.GetProp('_Name')] + list(map(lambda r: str(r[1]), val)))

    return 0


def submodule(name='__main__'):
    mdl = sys.modules[name]
    descs = get_descriptors_from_module(mdl)
    sys.exit(main(descs, prog=mdl.__spec__.name))


if __name__ == '__main__':
    sys.exit(main(all_descriptors(), prog='mordred'))
