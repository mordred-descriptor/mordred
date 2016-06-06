import six
import os
import sys
import click
import csv
from importlib import import_module
from . import __version__, all_descriptors, Calculator
from .error import MissingValue
from ._base import get_descriptors_from_module
from rdkit import Chem
from logging import getLogger
from tqdm import tqdm


logger = getLogger(__name__)


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
            logger.warn('smiles read failure: %s', name)
            continue

        mol.SetProp('_Name', name)
        yield mol


def sdf_parser(path):
    base = os.path.splitext(os.path.basename(path))[0]

    for i, mol in enumerate(Chem.SDMolSupplier(path, removeHs=False)):
        if mol is None:
            logger.warn('mol read failure: %s.%s', base, i)
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
        logger.warn('cannot detect file format: %s', path)
        r = ()

    for m in r:
        yield m


def callback_input(cxt, param, value):
    if len(value) == 0:
        if sys.stdin.isatty():
            click.echo(cxt.get_help())
            cxt.exit()
        else:
            return (sys.stdin,)

    return value


def callback_filetype(cxt, param, value):
    if value == 'auto':
        return auto_parser
    elif value == 'smi':
        return smiles_parser

    return sdf_parser


def callback_quiet(cxt, param, value):
    if cxt.params['output'].isatty():
        return True

    return value


@click.command(
    epilog='===== descriptors =====\n\n{}'.format(' '.join(all_descs)),
    context_settings={
        'help_option_names': ['-h', '--help']
    }
)
@click.version_option(
    __version__, '-v', '--version',
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
    default=sys.stdout, type=click.File('w') if six.PY3 else click.File('wb'),
    help='output csv file'
)
@click.option(
    'nproc', '-p', '--processes',
    default=None, type=click.INT,
    help='number of processes'
)
@click.option(
    '-q', '--quiet', callback=callback_quiet,
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
def main(input, parser, output, nproc, quiet, stream, descriptor, with3D):
    mols = (m for i in input for m in parser(i))

    if stream:
        N = None
    else:
        mols = list(mols)
        N = len(mols)

    # Descriptors
    calc = Calculator()
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
        writer = csv.writer(output)
        writer.writerow(['name'] + [str(d) for d in calc.descriptors])

        def warning(name, v, err_set):
            if not isinstance(v, MissingValue):
                return

            red = v.error.__class__, v.error.args
            if red in err_set:
                return

            tqdm.write('[{}] {}: {}'.format(v.header, name, v))
            err_set.add(red)

        def pretty(name, v, err_set):
            warning(name, v, err_set)

            if isinstance(v, MissingValue):
                return ''

            return str(v)

        for mol, val in calc.map(mols, nproc=nproc, nmols=N, quiet=quiet):
            err_set = set()

            if mol.HasProp('_Name'):
                name = mol.GetProp('_Name')
            else:
                name = Chem.MolToSmiles(mol)

            writer.writerow([name] + [pretty(name, v, err_set) for v in val])


if __name__ == '__main__':
    main(prog_name='{} -m {}'.format(os.path.basename(sys.executable), __package__))
