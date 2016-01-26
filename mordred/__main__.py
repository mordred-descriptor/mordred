import argparse
import csv
import math
import os
import sys

from rdkit import Chem

from ._base import Calculator, all_descriptors, get_descriptors_from_module


def smiles_parser(ifile):
    with ifile as f:
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


def dir_parser(d, fmt):
    for root, dirs, files in os.walk(d):
        for f in files:
            for m in file_parser(os.path.join(root, f), fmt):
                yield m


def file_parser(ifile, fmt):
    if fmt == 'smi':
        if isinstance(ifile, six.string_types):
            ifile = open(ifile)

        it = smiles_parser(ifile)
    elif fmt in ('sdf', 'mol'):
        it = sdf_parser(ifile)
    elif fmt == 'auto':
        if ifile == sys.stdin:
            it = file_parser(ifile, 'smi')
        elif os.path.isdir(ifile):
            it = dir_parser(ifile, fmt)
        else:
            ext = os.path.splitext(ifile)[1][1:]
            it = file_parser(ifile, ext)
    else:
        raise ValueError('unknown format: {!r}'.format(fmt))

    return it


def main(descs, prog=None):
    parser_options = dict()
    if prog is not None:
        parser_options['prog'] = '{} -m {}'.format(os.path.basename(sys.executable), prog)

    parser = argparse.ArgumentParser(**parser_options)

    parser.add_argument(
        '-i', '--input', metavar='PATH',
        type=str, default=sys.stdin,
        help='input file or directory(default: stdin)',
    )

    parser.add_argument(
        '-f', '--from', metavar='TYPE',
        default='auto', choices=['auto', 'smi', 'sdf', 'mol'],
        help='input filetype(one of %(choices)s, default: %(default)s)',
    )

    if sys.version_info >= (3, 0, 0):
        default_ofile = argparse.FileType('w')
    else:
        default_ofile = argparse.FileType('wb')

    parser.add_argument(
        '-o', '--output', metavar='PATH',
        type=default_ofile, default=sys.stdout,
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

        writer.writerow(['name'] + [str(d) for d in calc.descriptors])

        for mol, val in calc.map(mols, args.processes):
            def ppr(a):
                if isinstance(a, float) and math.isnan(a):
                    return ''
                else:
                    return str(a)

            writer.writerow([mol.GetProp('_Name')] + [ppr(v[1]) for v in val])

    return 0


def submodule(name='__main__'):
    mdl = sys.modules[name]
    descs = get_descriptors_from_module(mdl)

    prog = getattr(mdl, '__spac__', None)
    if prog is None:
        path, submdl = os.path.split(os.path.splitext(mdl.__file__)[0])
        pkg = os.path.basename(path)
        name = '{}.{}'.format(pkg, submdl)
    else:
        name = prog.name

    sys.exit(main(descs, prog=name))


if __name__ == '__main__':
    sys.exit(main(all_descriptors(), prog='mordred'))
