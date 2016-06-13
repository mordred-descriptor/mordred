from __future__ import print_function

import os
import sys
import logging
from importlib import import_module

import click
from rdkit import Chem

from . import Calculator, __version__, all_descriptors
from ._base import get_descriptors_from_module
from .error import Missing, MissingValueBase

all_descs = [
    '.'.join(m.__name__.split('.')[1:])
    for m in all_descriptors()
]


def smiles_parser(path):
    for line in open(path):
        line = line.strip().split()
        if len(line) == 1:
            smi = line[0]
            name = smi
        else:
            smi = line[0]
            name = ' '.join(line[1:])

        mol = Chem.MolFromSmiles(smi)

        if mol is None:
            logging.warning('smiles read failure: %s', name)
            continue

        mol.SetProp('_Name', name)
        yield mol


def sdf_parser(path):
    base = os.path.splitext(os.path.basename(path))[0]

    for i, mol in enumerate(Chem.SDMolSupplier(path, removeHs=False)):
        if mol is None:
            logging.warning('mol read failure: %s.%s', base, i)
            continue

        if mol.GetProp('_Name') == '':
            mol.SetProp(
                '_Name',
                '{}.{}'.format(base, i)
            )

        yield mol


def auto_parser(path):
    ext = os.path.splitext(path)[1]
    if ext == '.smi':
        r = smiles_parser(path)
    elif ext in ['.mol', '.sdf']:
        r = sdf_parser(path)
    else:
        logging.warning('cannot detect file format: %s', path)
        r = ()

    for m in r:
        yield m


def callback_input(cxt, param, value):
    if len(value) == 0:
        cxt.fail('INPUT file required')

    return value


def callback_filetype(cxt, param, value):
    if value == 'auto':
        return auto_parser
    elif value == 'smi':
        return smiles_parser

    return sdf_parser


@click.command(
    epilog='===== descriptors =====\n\n{}'.format(' '.join(all_descs)),
    context_settings={
        'help_option_names': ['-h', '--help']
    }
)
@click.version_option(
    __version__, '--version',
    prog_name='mordred'
)
@click.argument(
    'INPUT', nargs=-1,
    type=click.Path(exists=True),
    callback=callback_input
)
@click.option(
    'parser', '-t', '--type', default='auto',
    type=click.Choice(['auto', 'smi', 'mol', 'sdf']),
    help='input filetype',
    callback=callback_filetype
)
@click.option(
    '-o', '--output',
    default=sys.stdout, type=click.File('w'),
    help='output csv file'
)
@click.option(
    'nproc', '-p', '--processes',
    default=None, type=click.INT,
    help='number of processes'
)
@click.option(
    '-q', '--quiet',
    default=False, flag_value=True,
    help='hide progress bar'
)
@click.option(
    '-s', '--stream',
    default=False, flag_value=True,
    help='stream read'
)
@click.option(
    '-d', '--descriptor', multiple=True,
    type=click.Choice(all_descs),
    help='descriptors', metavar='DESC'
)
@click.option(
    'with3D', '-3', '--3D',
    default=False, flag_value=True,
    help='use 3D descriptors (require sdf or mol file)'
)
@click.option(
    '-v', '--verbosity',
    type=click.IntRange(0, 2), count=True,
    help='verbosity (0-2)'
)
def main(input, parser, output, nproc, quiet, stream, descriptor, with3D, verbosity):
    mols = (m for i in input for m in parser(i))

    if output.isatty():
        quiet = True

    if stream:
        N = None
    else:
        mols = list(mols)
        N = len(mols)

    # Descriptors
    calc = Calculator()

    if verbosity >= 2:
        calc._debug = True

    if len(descriptor) == 0:
        calc.register(all_descriptors(), exclude3D=not with3D)
    else:
        calc.register(
            (d
             for m in descriptor
             for d in get_descriptors_from_module(import_module('.' + m, __package__))
             ),
            exclude3D=not with3D
        )

    with output:
        write_row(output, ['name'] + [str(d) for d in calc.descriptors])

        def warning(name, v, err_set):
            if not isinstance(v, MissingValueBase):
                return

            if verbosity == 0 and isinstance(v, Missing):
                return

            red = v.error.__class__, v.error.args
            if red in err_set:
                return

            calc.echo('[{}] {}: {}'.format(v.header, name, v), file=sys.stderr)
            err_set.add(red)

        def pretty(name, v, err_set):
            warning(name, v, err_set)

            if isinstance(v, MissingValueBase):
                return ''

            return str(v)

        for mol, val in calc.map(mols, nproc=nproc, nmols=N, quiet=quiet):
            err_set = set()

            if mol.HasProp('_Name'):
                name = mol.GetProp('_Name')
            else:
                name = Chem.MolToSmiles(mol)

            write_row(output, [name] + [pretty(name, v, err_set) for v in val])


def write_row(file, data):
    file.write(
        ','.join(
            str(v).replace('"', '""').replace('\n', '').replace('\r', '')
            for v in data
        )
    )
    file.write('\n')


if __name__ == '__main__':
    main(prog_name='{} -m {}'.format(os.path.basename(sys.executable), __package__))
