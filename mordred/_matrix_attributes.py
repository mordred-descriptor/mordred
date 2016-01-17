from ._base import Descriptor
import numpy as np
from collections import namedtuple
from six import string_types


Eig = namedtuple('eigen', 'val vec min max')

methods = []


def method(cls):
    methods.append(cls)
    return cls


class common(Descriptor):
    descriptor_keys = 'matrix', 'explicit_hydrogens', 'gasteiger_charges', 'kekulize'

    def __init__(self, matrix, explicit_hydrogens=True, gasteiger_charges=False, kekulize=False):
        self.matrix = matrix
        self.explicit_hydrogens = explicit_hydrogens
        self.gasteiger_charges = gasteiger_charges
        self.kelulize = kekulize

    @property
    def _key_args(self):
        return [
            self.matrix,
            self.explicit_hydrogens,
            self.gasteiger_charges,
            self.kekulize
        ]

    @property
    def _eig(self):
        return Eigen(*self._key_args)

    @property
    def _SpMax(self):
        return SpMax(*self._key_args)

    @property
    def _SpMean(self):
        return SpMean(*self._key_args)

    @property
    def _SpAD(self):
        return SpAD(*self._key_args)

    @property
    def _VE1(self):
        return VE1(*self._key_args)

    @property
    def _VR1(self):
        return VR1(*self._key_args)


class Eigen(common):
    @property
    def dependencies(self):
        return dict(matrix=self.matrix)

    def calculate(self, mol, matrix):
        w, v = np.linalg.eig(matrix)

        if np.iscomplexobj(w):
            w = w.real

        if np.iscomplexobj(v):
            v = v.real

        i_min = np.argmin(w)
        i_max = np.argmax(w)

        return Eig(w, v, i_min, i_max)


@method
class SpAbs(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        return np.abs(eig.val).sum()


@method
class SpMax(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        return eig.val[eig.max]


@method
class SpDiam(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig, SpMax=self._SpMax)

    def calculate(self, mol, SpMax, eig):
        return SpMax - eig.val[eig.min]


class SpMean(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        return np.mean(eig.val)


@method
class SpAD(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig, SpMean=self._SpMean)

    def calculate(self, mol, eig, SpMean):
        return np.abs(eig.val - SpMean).sum()


@method
class SpMAD(common):
    @property
    def dependencies(self):
        return dict(SpAD=self._SpAD)

    def calculate(self, mol, SpAD):
        return SpAD / mol.GetNumAtoms()


@method
class LogEE(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        # log sum exp: https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp
        a = np.maximum(eig.val[eig.max], 0)
        sx = np.exp(eig.val - a).sum() + np.exp(-a)
        return a + np.log(sx)


@method
class SM1(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        return eig.val.sum()


@method
class VE1(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        return np.abs(eig.vec[:, eig.max]).sum()


@method
class VE2(common):
    @property
    def dependencies(self):
        return dict(VE1=self._VE1)

    def calculate(self, mol, VE1):
        return VE1 / mol.GetNumAtoms()


@method
class VE3(common):
    @property
    def dependencies(self):
        return dict(VE1=self._VE1)

    def calculate(self, mol, VE1):
        if VE1 == 0:
            return np.nan
        else:
            return np.log(0.1 * mol.GetNumAtoms() * VE1)


@method
class VR1(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        s = 0

        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            s += (eig.vec[i, eig.max] * eig.vec[j, eig.max]) ** -0.5

        return s


@method
class VR2(common):
    @property
    def dependencies(self):
        return dict(VR1=self._VR1)

    def calculate(self, mol, VR1):
        return VR1 / mol.GetNumAtoms()


@method
class VR3(common):
    @property
    def dependencies(self):
        return dict(VR1=self._VR1)

    def calculate(self, mol, VR1):
        if VR1 == 0:
            return np.nan
        else:
            return np.log(0.1 * mol.GetNumAtoms() * VR1)


method_dict = {m.__name__: m for m in methods}


def get_method(n):
    if isinstance(n, string_types):
        return method_dict[n]

    return n
