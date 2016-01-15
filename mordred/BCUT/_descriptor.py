from .._base import Descriptor
from .. import _atomic_property
import numpy as np


class BurdenMatrixDescriptor(Descriptor):
    explicit_hydrogens = False


class burden(BurdenMatrixDescriptor):
    def calculate(self, mol):
        N = mol.GetNumAtoms()

        mat = 0.001 * np.ones((N, N))

        for bond in mol.GetBonds():
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()
            i = a.GetIdx()
            j = b.GetIdx()

            try:
                w = bond.GetBondTypeAsDouble() / 10.0
            except RuntimeError:
                w = 1.0

            if a.GetDegree() == 1 or b.GetDegree() == 1:
                w += 0.01

            mat[i, j] = w
            mat[j, i] = w

        return mat


class burden_eigen_values(BurdenMatrixDescriptor):
    @property
    def descriptor_key(self):
        return self.make_key(self.prop, self.gasteiger_charges)

    def __init__(self, prop, gasteiger_charges):
        self.prop = prop
        self.gasteiger_charges = gasteiger_charges

    @property
    def dependencies(self):
        return dict(burden=burden.make_key())

    def calculate(self, mol, burden):
        bmat = burden.copy()
        np.fill_diagonal(bmat, [self.prop(a) for a in mol.GetAtoms()])
        return np.linalg.eig(bmat)[0]


class BCUT(BurdenMatrixDescriptor):
    r'''
    BCUT descriptor

    Parameters:
        prop(str, function): atomic property
        from_high(bool): n-th eigen value from high
        nth(int): n-th eigen value

    Returns:
        float: result
    '''

    @classmethod
    def preset(cls):
        return (cls(a, 1, h) for a in _atomic_property.get_properties(istate=True) for h in [False, True])

    @property
    def descriptor_name(self):
        hl = 'h' if self.from_high else 'l'
        return 'BCUT{}-{}{}'.format(self.prop_name, self.nth, hl)

    @property
    def descriptor_key(self):
        return self.make_key(self.prop, self.nth, self.from_high)

    def __init__(self, prop='m', nth=1, from_high=True):
        if prop == 'c':
            self.prop_name = 'c'
            self.prop = _atomic_property.get_charge_implicitHs
            self.gasteiger_charges = True
        else:
            self.prop_name, self.prop = _atomic_property.getter(prop)
            self.gasteiger_charges = False

        self.nth = nth
        self.from_high = from_high

    @property
    def dependencies(self):
        return dict(bev=burden_eigen_values.make_key(self.prop, self.gasteiger_charges))

    def calculate(self, mol, bev):
        nth = self.nth - 1
        if 0 <= nth < len(bev):
            if self.from_high:
                nth = len(bev) - nth - 1
            return np.sort(bev)[nth]
        else:
            return np.nan


_descriptors = [BCUT]
__all__ = [d.__name__ for d in _descriptors]
