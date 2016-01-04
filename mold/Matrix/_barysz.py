from .._base import Descriptor
from .. import _atomic_property
from rdkit import Chem
from scipy.sparse.csgraph import shortest_path
import numpy as np
from ._matrix_attributes import methods


class BaryszMatrixDescriptor(Descriptor):
    explicit_hydrogens = False


class barysz(BaryszMatrixDescriptor):
    _carbon = Chem.Atom(6)

    @property
    def descriptor_key(self):
        return self.make_key(self.prop)

    def __init__(self, prop):
        self.prop = prop

    def calculate(self, mol):
        N = mol.GetNumAtoms()

        C = self.prop(self._carbon)

        dmat = np.zeros((N, N))

        for bond in mol.GetBonds():
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()

            i = a.GetIdx()
            j = b.GetIdx()

            pi = bond.GetBondTypeAsDouble()

            w = float(C * C) / float(self.prop(a) * self.prop(b) * pi)
            dmat[i, j] = w
            dmat[j, i] = w

        sp = shortest_path(dmat)
        np.fill_diagonal(sp, [1. - float(C) / self.prop(a) for a in mol.GetAtoms()])
        return sp


method_dict = {m.__name__: m for m in methods}


class BaryszMatrix(BaryszMatrixDescriptor):
    descriptor_defaults = [(p, m.__name__) for p in 'Zmvepi' for m in methods]

    @property
    def descriptor_key(self):
        return self.make_key(self.prop, self.method)

    @property
    def dependencies(self):
        return dict(result=method_dict[self.method].make_key(
            barysz.make_key(self.prop),
            self.explicit_hydrogens,
            self.gasteiger_charges,
            self.kekulize,
        ))

    @property
    def descriptor_name(self):
        return '{}_Dz{}'.format(self.method, self.prop_name)

    def __init__(self, prop, method):
        self.prop_name, self.prop = _atomic_property.getter(prop)
        self.method = method

    def calculate(self, mol, result):
        return result
