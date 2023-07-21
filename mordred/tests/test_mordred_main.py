import os
import sys
import shutil
import tempfile
from contextlib import contextmanager

from rdkit.Chem import AllChem as Chem


from .. import Calculator, descriptors
from ..__main__ import main as mordred

Nd2D = len(Calculator(descriptors, ignore_3D=True).descriptors)
Nd3D = len(Calculator(descriptors, ignore_3D=False).descriptors)

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version

__version__ = version("mordredcommunity")


def in_(a, s):
    assert a in s, "{!r} not in {!r}".format(a, s)


@contextmanager
def isolate():
    cwd = os.getcwd()
    try:
        work = tempfile.mkdtemp()
        os.chdir(work)
        yield
    finally:
        os.chdir(cwd)
        shutil.rmtree(work)


class Result(object):
    def __iter__(self):
        return iter([self.stdout, self.stderr, self.exitcode])

    def __repr__(self):
        return "{!r} {!r} {!r}".format(self.stdout, self.stderr, self.exitcode)


def command(cmd, *args):
    r = Result()

    orig = sys.stdout, sys.stderr

    ofd, opath = tempfile.mkstemp()
    efd, epath = tempfile.mkstemp()

    try:
        with os.fdopen(ofd, "w") as stdout, os.fdopen(efd, "w") as stderr:
            sys.stdout = stdout
            sys.stderr = stderr

            try:
                cmd(args=args)
                r.exitcode = 0
            except SystemExit as e:
                r.exitcode = e.args[0]

    finally:
        sys.stdout, sys.stderr = orig

        with open(opath) as stdout, open(epath) as stderr:
            r.stdout = stdout.read()
            r.stderr = stderr.read()

        os.remove(opath)
        os.remove(epath)

    return r


def test_no_args():
    stdout, stderr, exitcode = command(mordred)
    assert exitcode == 2
    assert stdout == ""
    in_("usage:", stderr)
    # python3 or python2
    assert (
        "the following arguments are required: INPUT" in stderr
        or "too few arguments" in stderr
    )


def test_help():
    stdout, stderr, exitcode = command(mordred, "-h")
    assert exitcode == 0
    assert stderr == ""
    in_("usage:", stdout)
    in_("descriptors:", stdout)


def test_version():
    stdout, stderr, exitcode = command(mordred, "--version")
    assert exitcode == 0

    vstr = "mordred-{}\n".format(__version__)

    assert stdout == vstr


def test_missing_file():
    with isolate():
        stdout, stderr, exitcode = command(mordred, "missing.smi")
        assert exitcode == 2
        assert stdout == ""
        in_("usage:", stderr)
        in_("invalid PathType value", stderr)


def number_of_field(output):
    return len(list(filter(lambda c: c == ",", output.split("\n")[0])))


def test_smi():
    with isolate():
        with open("input.smi", "w") as f:
            f.write("c1ccccc1 Benzene\n")

        stdout, stderr, exitcode = command(mordred, "input.smi", "-q", "-o", "-")

        assert exitcode == 0
        assert number_of_field(stdout) == Nd2D
        assert stdout.split("\n")[1].split(",")[0] == "Benzene"
        assert stderr == ""


def test_smi_without_name():
    with isolate():
        with open("input.smi", "w") as f:
            f.write("c1ccccc1\n")

        stdout, stderr, exitcode = command(mordred, "input.smi", "-q", "-o", "-")

        assert exitcode == 0
        assert number_of_field(stdout) == Nd2D
        assert stdout.split("\n")[1].split(",")[0] == "c1ccccc1"
        assert stderr == ""


def test_sdf():
    with isolate():
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol.SetProp("_Name", "Benzene")
        Chem.MolToMolFile(mol, "input.sdf")

        stdout, stderr, exitcode = command(
            mordred, "input.sdf", "-q", "-o", "output.csv"
        )

        assert exitcode == 0
        assert stdout == ""
        assert stderr == ""
        output = open("output.csv").read()
        assert number_of_field(output) == Nd2D
        assert output.split("\n")[1].split(",")[0] == "Benzene"


def test_sdf_3D():
    with isolate():
        mol = Chem.MolFromSmiles("c1ccccc1")
        Chem.EmbedMolecule(mol)
        mol.SetProp("_Name", "Benzene")
        Chem.MolToMolFile(mol, "input.sdf")

        stdout, stderr, exitcode = command(
            mordred, "input.sdf", "-q", "-o", "output.csv", "-3"
        )

        assert exitcode == 0
        assert stdout == ""
        assert stderr == ""
        output = open("output.csv").read()
        assert number_of_field(output) == Nd3D
        assert output.split("\n")[1].split(",")[0] == "Benzene"


def test_verbose():
    with isolate():
        with open("input.smi", "w") as f:
            f.write("c1ccccc1 Benzene\n")

        r = command(mordred, "input.smi", "-o", "-", "-q", "-vv")

        assert r.exitcode == 0
        assert number_of_field(r.stdout) == Nd2D
        assert r.stdout.split("\n")[1].split(",")[0] == "Benzene"
        assert "[Missing]" in r.stderr
