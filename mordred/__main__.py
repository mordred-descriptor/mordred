from __future__ import print_function

import os
import sys
import logging
import argparse
from importlib import import_module
from multiprocessing import freeze_support

from rdkit import Chem

from . import Calculator, descriptors
from ._base import get_descriptors_in_module
from ._util import PathType, module_prog
from .error import Missing, MissingValueBase

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version

__version__ = version("mordredcommunity")


def smiles_parser(path):
    with open(path) as file:
        for line in file:
            line = line.strip().split()
            if len(line) == 1:
                smi = line[0]
                name = smi
            else:
                smi = line[0]
                name = " ".join(line[1:])

            mol = Chem.MolFromSmiles(smi)

            if mol is None:
                logging.warning("smiles read failure: %s", name)
                continue

            mol.SetProp("_Name", name)
            yield mol


def sdf_parser(path):
    base = os.path.splitext(os.path.basename(path))[0]

    for i, mol in enumerate(Chem.SDMolSupplier(path, removeHs=False)):
        if mol is None:
            logging.warning("mol read failure: %s.%s", base, i)
            continue

        if mol.GetProp("_Name") == "":
            mol.SetProp("_Name", "{}.{}".format(base, i))

        yield mol


def auto_parser(path):
    ext = os.path.splitext(path)[1]
    if ext == ".smi":
        r = smiles_parser(path)
    elif ext in [".mol", ".sdf"]:
        r = sdf_parser(path)
    else:
        logging.warning("cannot detect file format: %s", path)
        r = ()

    for m in r:
        yield m


class ParserAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super(ParserAction, self).__init__(option_strings, dest, **kwargs)
        self.default = self.to_parser("auto")
        self.choices = ["auto", "sdf", "mol", "smi"]

    def to_parser(self, value):
        if value == "auto":
            return auto_parser
        elif value == "smi":
            return smiles_parser
        elif value in ["sdf", "mol"]:
            return sdf_parser

        raise ValueError("invalid parser: {}".format(value))

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, self.to_parser(values))


def make_parser():
    parser = argparse.ArgumentParser(
        prog=module_prog(__package__),
        epilog="descriptors: {}".format(" ".join(descriptors.__all__)),
    )
    parser.add_argument(
        "--version",
        action="version",
        help="input molecular file",
        version="{}-{}".format(__package__, __version__),
    )
    parser.add_argument("input", type=PathType, nargs="+", metavar="INPUT")
    parser.add_argument(
        "-t", "--type", action=ParserAction, help="input filetype (default: auto)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        type=argparse.FileType("w"),
        help="output file path (default: stdout)",
    )
    parser.add_argument(
        "-p",
        "--processes",
        default=None,
        type=int,
        help="number of processes (default: number of logical processors)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="hide progress bar")
    parser.add_argument("-s", "--stream", action="store_true", help="stream read")
    parser.add_argument(
        "-d",
        "--descriptor",
        default=[],
        choices=descriptors.__all__,
        action="append",
        help="descriptors to calculate (default: all)",
        metavar="DESC",
    )
    parser.add_argument(
        "-3",
        "--3D",
        action="store_true",
        dest="with3D",
        help="use 3D descriptors (require sdf or mol file)",
    )
    parser.add_argument(
        "-v", "--verbosity", action="count", default=0, help="verbosity"
    )

    return parser


def main_process(
    input, parser, output, nproc, quiet, stream, descriptor, with3D, verbosity
):
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
        calc.register(descriptors, ignore_3D=not with3D)
    else:
        calc.register(
            (
                d
                for m in descriptor
                for d in get_descriptors_in_module(
                    import_module("." + m, __package__), False
                )
            ),
            ignore_3D=not with3D,
        )

    with output:
        write_row(output, ["name"] + [str(d) for d in calc.descriptors])

        def warning(name, v, err_set):
            if not isinstance(v, MissingValueBase):
                return

            if verbosity == 0 and isinstance(v, Missing):
                return

            red = v.error.__class__, v.error.args
            if red in err_set:
                return

            calc.echo("[{}] {}: {}".format(v.header, name, v), file=sys.stderr)
            err_set.add(red)

        def pretty(name, v, err_set):
            warning(name, v, err_set)

            if isinstance(v, MissingValueBase):
                return ""

            if isinstance(v, bool):
                return int(v)

            return str(v)

        for result in calc.map(mols, nproc=nproc, nmols=N, quiet=quiet):
            err_set = set()

            if result.mol.HasProp("_Name"):
                name = result.mol.GetProp("_Name")
            else:
                name = Chem.MolToSmiles(result.mol)

            write_row(output, [name] + [pretty(name, v, err_set) for v in result])


def write_row(file, data):
    file.write(
        ",".join(
            str(v).replace('"', '""').replace("\n", "").replace("\r", "") for v in data
        )
    )
    file.write("\n")


def main(args=None):
    parser = make_parser()
    p = parser.parse_args(args)
    return main_process(
        input=p.input,
        parser=p.type,
        output=p.output,
        nproc=p.processes,
        quiet=p.quiet,
        stream=p.stream,
        descriptor=p.descriptor,
        with3D=p.with3D,
        verbosity=p.verbosity,
    )


if __name__ == "__main__":
    freeze_support()
    main()
