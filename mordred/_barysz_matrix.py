from ._base import Descriptor
from . import _atomic_property
from rdkit import Chem
from scipy.sparse.csgraph import shortest_path
import numpy as np
from ._matrix_attributes import methods, get_method


class BaryszMatrixBase(Descriptor):
    explicit_hydrogens = False

    @property
    def gasteiger_charges(self):
        return getattr(self.prop, 'gasteiger_charges', False)


_carbon = Chem.Atom(6)


class barysz(BaryszMatrixBase):
    descriptor_keys = 'prop',

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


class BaryszMatrix(BaryszMatrixBase):
    r'''
    barysz matrix descriptor

    Parameters:
        prop(str, function): atomic property
        type(str): matrix aggregateing method

    Returns:
        float: result
    '''

    @classmethod
    def preset(cls):
        return (cls(p, m) for p in _atomic_property.get_properties() for m in methods)

    def __str__(self):
        return '{}_Dz{}'.format(self.type.__name__, self.prop_name)

    descriptor_keys = 'prop', 'type'

    def __init__(self, prop='Z', type='SpMax'):
        self.prop_name, self.prop = _atomic_property.getter(prop, self.explicit_hydrogens)
        self.type = get_method(type)

    @property
    def dependencies(self):
        return dict(
            result=self.type(
                barysz(self.prop),
                self.explicit_hydrogens,
                self.gasteiger_charges,
                self.kekulize,
            )
        )

    def calculate(self, mol, result):
        return result
