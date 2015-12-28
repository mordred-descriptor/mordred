from .._base import Descriptor
from .. import _atomic_property
import numpy as np


class BurdenMatrixDescriptor(Descriptor):
    explicitHydrogens = False


class burden(BurdenMatrixDescriptor):
    @property
    def descriptor_key(self):
        return self.make_key()

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
                continue

            if a.GetDegree() == 1 or b.GetDegree() == 1:
                w += 0.01

            mat[i, j] = w
            mat[j, i] = w

        return mat


class burden_eigen_values(BurdenMatrixDescriptor):
    @property
    def gasteigerCharges(self):
        return self.charge

    @property
    def descriptor_key(self):
        return self.make_key(self.prop, self.charge)

    def __init__(self, prop, charge):
        self.prop = prop
        self.charge = charge

    @property
    def dependencies(self):
        return dict(burden=burden.make_key())

    def calculate(self, mol, burden):
        bmat = burden.copy()
        np.fill_diagonal(bmat, [self.prop(a) for a in mol.GetAtoms()])
        return np.linalg.eig(bmat)[0]


class BCUT(BurdenMatrixDescriptor):
    descriptor_defaults = [('m', 1, False), ('m', 1, True),
                           ('c', 1, False), ('c', 1, True)]

    @property
    def gasteigerCharges(self):
        return self.charge

    @property
    def descriptor_name(self):
        hl = 'h' if self.from_high else 'l'
        return 'BCUT{}-{}{}'.format(self.prop_name, self.nth, hl)

    @property
    def descriptor_key(self):
        return self.make_key(self.prop, self.nth, self.from_high)

    def __init__(self, prop, nth, from_high):
        if prop == 'c':
            self.prop_name = 'c'
            self.prop = _atomic_property.get_charge_implicitHs
            self.charge = True
        else:
            self.prop_name, self.prop = _atomic_property.getter(prop)
            self.charge = False

        self.nth = nth
        self.from_high = from_high

    @property
    def dependencies(self):
        return dict(bev=burden_eigen_values.make_key(self.prop, self.charge))

    def calculate(self, mol, bev):
        nth = self.nth - 1
        if 0 <= nth < len(bev):
            if self.from_high:
                nth = len(bev) - nth - 1
            return np.sort(bev)[nth]
        else:
            return np.nan
