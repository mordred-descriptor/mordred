from __future__ import print_function

import os
import sys

import numpy as np


def parse_enum(enum, v):
    if isinstance(v, enum):
        return v
    else:
        return enum[v]


def atoms_to_numpy(f, mol, dtype="float"):
    return np.fromiter((f(a) for a in mol.GetAtoms()), dtype, mol.GetNumAtoms())


def conformer_to_numpy(conf):
    return np.array([list(conf.GetAtomPosition(i)) for i in range(conf.GetNumAtoms())])


class Capture(object):
    def __init__(self, target="stderr"):
        self.target = target
        self.orig = getattr(sys, target)
        self.result = []

    def write(self, text):
        self.result.append(text)

    def flush(self):
        pass

    def __enter__(self):
        setattr(sys, self.target, self)
        return self

    def __exit__(self, *args):
        setattr(sys, self.target, self.orig)


class DummyBar(object):
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args, **kws):
        pass

    def update(self, *args, **kws):
        pass

    @classmethod
    def write(cls, s, file=sys.stdout, end="\n"):
        print(s, file=file, end=end)  # noqa: T003


class NotebookWrapper(object):
    def __init__(self, **kwargs):
        from tqdm import tqdm_notebook

        self.bar = tqdm_notebook(**kwargs)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def update(self, *args, **kwargs):
        self.bar.update(*args, **kwargs)

    def write(self, *args, **kwargs):
        self.bar.update(*args, **kwargs)


def PathType(string):
    if not os.path.isfile(string):
        raise ValueError("file not exists: {}".format(string))

    return string


def module_prog(pkg):
    return "{} -m {}".format(os.path.basename(sys.executable), pkg)


def to_ordinal(n):
    r"""Int to ordinal string.

    >>> to_ordinal(1)
    'first'
    >>> to_ordinal(2)
    'second'
    >>> to_ordinal(3)
    'third'
    >>> to_ordinal(4)
    '4-th'
    >>> to_ordinal(104)
    '104-th'

    """
    if n == 1:
        return "first"
    elif n == 2:
        return "second"
    elif n == 3:
        return "third"
    else:
        return "{}-th".format(n)
