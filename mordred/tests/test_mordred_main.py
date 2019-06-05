import os
import sys
import shutil
import tempfile
from contextlib import contextmanager

from rdkit.Chem import AllChem as Chem

from nose.tools import eq_

from .. import Calculator, __version__, descriptors
from ..__main__ import main as mordred

Nd2D = len(Calculator(descriptors, ignore_3D=True).descriptors)
Nd3D = len(Calculator(descriptors, ignore_3D=False).descriptors)


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
    eq_(exitcode, 2)
    eq_(stdout, "")
    in_("usage:", stderr)
    # python3 or python2
    assert (
        "the following arguments are required: INPUT" in stderr
        or "too few arguments" in stderr
    )


def test_help():
    stdout, stderr, exitcode = command(mordred, "-h")
    eq_(exitcode, 0)
    eq_(stderr, "")
    in_("usage:", stdout)
    in_("descriptors:", stdout)


def test_version():
    stdout, stderr, exitcode = command(mordred, "--version")
    eq_(exitcode, 0)

    vstr = "mordred-{}\n".format(__version__)

    if stderr == "":  # python 3
        eq_(stdout, vstr)
    else:  # python2
        eq_(stderr, vstr)


def test_missing_file():
    with isolate():
        stdout, stderr, exitcode = command(mordred, "missing.smi")
        eq_(exitcode, 2)
        eq_(stdout, "")
        in_("usage:", stderr)
        in_("invalid PathType value", stderr)


def number_of_field(output):
    return len(list(filter(lambda c: c == ",", output.split("\n")[0])))


def test_smi():
    with isolate():
        with open("input.smi", "w") as f:
            f.write("c1ccccc1 Benzene\n")

        stdout, stderr, exitcode = command(mordred, "input.smi", "-q", "-o", "-")

        eq_(exitcode, 0)
        eq_(number_of_field(stdout), Nd2D)
        eq_(stdout.split("\n")[1].split(",")[0], "Benzene")
        eq_(stderr, "")


def test_smi_without_name():
    with isolate():
        with open("input.smi", "w") as f:
            f.write("c1ccccc1\n")

        stdout, stderr, exitcode = command(mordred, "input.smi", "-q", "-o", "-")

        eq_(exitcode, 0)
        eq_(number_of_field(stdout), Nd2D)
        eq_(stdout.split("\n")[1].split(",")[0], "c1ccccc1")
        eq_(stderr, "")


def test_sdf():
    with isolate():
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol.SetProp("_Name", "Benzene")
        Chem.MolToMolFile(mol, "input.sdf")

        stdout, stderr, exitcode = command(
            mordred, "input.sdf", "-q", "-o", "output.csv"
        )

        eq_(exitcode, 0)
        eq_(stdout, "")
        eq_(stderr, "")
        output = open("output.csv").read()
        eq_(number_of_field(output), Nd2D)
        eq_(output.split("\n")[1].split(",")[0], "Benzene")


def test_sdf_3D():
    with isolate():
        mol = Chem.MolFromSmiles("c1ccccc1")
        Chem.EmbedMolecule(mol)
        mol.SetProp("_Name", "Benzene")
        Chem.MolToMolFile(mol, "input.sdf")

        stdout, stderr, exitcode = command(
            mordred, "input.sdf", "-q", "-o", "output.csv", "-3"
        )

        eq_(exitcode, 0)
        eq_(stdout, "")
        eq_(stderr, "")
        output = open("output.csv").read()
        eq_(number_of_field(output), Nd3D)
        eq_(output.split("\n")[1].split(",")[0], "Benzene")


def test_verbose():
    with isolate():
        with open("input.smi", "w") as f:
            f.write("c1ccccc1 Benzene\n")

        r = command(mordred, "input.smi", "-o", "-", "-q", "-vv")

        eq_(r.exitcode, 0)
        eq_(number_of_field(r.stdout), Nd2D)
        eq_(r.stdout.split("\n")[1].split(",")[0], "Benzene")
        in_("[Missing]", r.stderr)
