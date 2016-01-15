from ._base import Descriptor
from . import _atomic_property
from rdkit import Chem
from scipy.sparse.csgraph import shortest_path
import numpy as np
from ._matrix_attributes import methods, get_method


class BaryszMatrixDescriptor(Descriptor):
    explicit_hydrogens = False


_carbon = Chem.Atom(6)


class barysz(BaryszMatrixDescriptor):
    @property
    def descriptor_key(self):
        return self.make_key(self.prop)

    def __init__(self, prop):
        self.prop = prop

    def calculate(self, mol):
        N = mol.GetNumAtoms()

        C = self.prop(_carbon)

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


class BaryszMatrix(BaryszMatrixDescriptor):
    r'''
    barysz matrix descriptor

    Parameters:
        prop(str, function): atomic property
        method(str): matrix aggregate method

    Returns:
        float: result
    '''

    @classmethod
    def preset(cls):
        return (cls(p, m) for p in _atomic_property.get_properties() for m in methods)

    @property
    def descriptor_key(self):
        return self.make_key(self.prop, self.method)

    @property
    def dependencies(self):
        return dict(result=self.method.make_key(
            barysz.make_key(self.prop),
            self.explicit_hydrogens,
            self.gasteiger_charges,
            self.kekulize,
        ))

    @property
    def descriptor_name(self):
        return '{}_Dz{}'.format(self.method.__name__, self.prop_name)

    def __init__(self, prop='Z', method='SpMax'):
        self.prop_name, self.prop = _atomic_property.getter(prop)
        self.method = get_method(method)

    def calculate(self, mol, result):
        return result
