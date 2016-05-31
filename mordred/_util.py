import numpy as np


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
