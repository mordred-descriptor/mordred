import sys
import numpy as np
from tqdm import tqdm


def parse_enum(enum, v):
    if isinstance(v, enum):
        return v
    else:
        return enum[v]


def atoms_to_numpy(f, mol, dtype='float'):
    return np.fromiter(
        (f(a) for a in mol.GetAtoms()),
        dtype, mol.GetNumAtoms()
    )


def conformer_to_numpy(conf):
    return np.array(
        [list(conf.GetAtomPosition(i)) for i in range(conf.GetNumAtoms())]
    )


class Capture(object):
    def __init__(self, target='stderr'):
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
    def __init__(self, logger):
        self.logger = logger

    def __enter__(self):
        return self

    def __exit__(self, *args, **kws):
        pass

    def update(self, *args, **kws):
        pass

    def write(self, text, *args, **kwargs):
        self.logger.warn(text)


def get_bar(quiet, logger, total, **kwargs):
    if quiet:
        return DummyBar(logger)
    else:
        args = dict(
            dynamic_ncols=True,
            leave=True,
        )
        args.update(kwargs)
        args['total'] = total
        return tqdm(**args)
