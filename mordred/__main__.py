import argparse
import csv
import math
import os
import sys

from rdkit import Chem

import six

import tqdm

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
                sys.stderr.flush()
                continue

            mol.SetProp('_Name', name)
            yield mol


def sdf_parser(path):
    for mol in Chem.SDMolSupplier(path, removeHs=False):
        yield mol


def mol_parser(path):
    mol = Chem.MolFromMolFile(path)
    if mol is None:
        return ()

    if mol.GetProp('_Name') == '':
        mol.SetProp(
            '_Name',
            os.path.splitext(os.path.basename(path))[0],
        )

    return [mol]


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
    elif fmt in 'sdf':
        it = sdf_parser(ifile)
    elif fmt == 'mol':
        it = mol_parser(ifile)
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

    return (mol for mol in it if mol is not None)


class DummyBar(object):
    def __init__(self, *args, **kws):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args, **kws):
        pass

    def update(self, *args, **kws):
        pass


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

    parser.add_argument(
        '-o', '--output', metavar='PATH',
        type=str, default=None,
        help='output csv file(default: stdout)'
    )

    parser.add_argument(
        '-p', '--processes', metavar='N',
        default=None, type=int,
        help='number of processes to use(default: number of threads)',
    )

    parser.add_argument(
        '-q', '--quiet',
        default=False, action='store_true',
        help="hide progress bar",
    )

    parser.add_argument(
        '-s', '--stream',
        default=False, action='store_true',
        help='stream read',
    )

    args = parser.parse_args()

    if args.input == sys.stdin and args.input.isatty():
        sys.exit(parser.print_help())

    if args.output is None:
        args.output = sys.stdout
    elif six.PY3:
        args.output = open(args.output, 'w', newline='')
    else:
        args.output = open(args.output, 'wb')

    if args.output.isatty():
        args.quiet = True

    mols = file_parser(args.input, getattr(args, 'from'))

    if not args.stream:
        mols = list(mols)
        N = len(mols)
    else:
        N = None

    calc = Calculator(descs)

    with args.output as output:
        progress_args = {'dynamic_ncols': True, 'leave': True}
        if args.quiet:
            Progress = DummyBar
        elif N is None:
            Progress = tqdm.tqdm
        else:
            Progress = tqdm.tqdm
            progress_args.update({'total': N})

        writer = csv.writer(output)

        writer.writerow(['name'] + [str(d) for d in calc.descriptors])

        with Progress(**progress_args) as bar:
            def callback(r):
                bar.update()

            for mol, val in calc.map(
                    mols, args.processes, error_mode='log', callback=callback):

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
