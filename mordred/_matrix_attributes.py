from collections import namedtuple

import numpy as np

from six import string_types

from ._base import Descriptor
from .exception import MordredValueError

Eig = namedtuple('eigen', 'val vec min max')

methods = []


def method(cls):
    methods.append(cls)
    cls.as_argument = cls.__name__
    return cls


class Common(Descriptor):
    require_connected = True

    def as_key(self):
        return self.__class__, (
            self.matrix,
            self.explicit_hydrogens,
            self.kekulize,
        )

    def __init__(self, matrix, explicit_hydrogens, kekulize):
        self.matrix = matrix
        self.explicit_hydrogens = explicit_hydrogens
        self.kelulize = kekulize

    @property
    def _key_args(self):
        return (
            self.matrix,
            self.explicit_hydrogens,
            self.kekulize
        )

    def dependencies(self):
        return {'eig': Eigen(*self._key_args)}

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


class Eigen(Common):
    def dependencies(self):
        return {'matrix': self.matrix}

    def calculate(self, matrix):
        if matrix is None:
            raise MordredValueError('Eigen: matrix is None')

        w, v = np.linalg.eig(matrix)

        if np.iscomplexobj(w):
            w = w.real

        if np.iscomplexobj(v):
            v = v.real

        i_min = np.argmin(w)
        i_max = np.argmax(w)

        return Eig(w, v, i_min, i_max)


@method
class SpAbs(Common):
    def calculate(self, eig):
        return np.abs(eig.val).sum()


@method
class SpMax(Common):
    def calculate(self, eig):
        return eig.val[eig.max]


@method
class SpDiam(Common):
    def dependencies(self):
        return {
            'eig': self._eig,
            'SpMax': self._SpMax,
        }

    def calculate(self, SpMax, eig):
        return SpMax - eig.val[eig.min]


class SpMean(Common):
    def calculate(self, eig):
        return np.mean(eig.val)


@method
class SpAD(Common):
    def dependencies(self):
        return {
            'eig': self._eig,
            'SpMean': self._SpMean,
        }

    def calculate(self, eig, SpMean):
        return np.abs(eig.val - SpMean).sum()


@method
class SpMAD(Common):
    def dependencies(self):
        return {'SpAD': self._SpAD}

    def calculate(self, SpAD):
        return SpAD / self.mol.GetNumAtoms()


@method
class LogEE(Common):
    def calculate(self, eig):
        # log sum exp: https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp
        a = np.maximum(eig.val[eig.max], 0)
        sx = np.exp(eig.val - a).sum() + np.exp(-a)
        return a + np.log(sx)


@method
class SM1(Common):
    def calculate(self, eig):
        return eig.val.sum()


@method
class VE1(Common):
    def calculate(self, eig):
        return np.abs(eig.vec[:, eig.max]).sum()


@method
class VE2(Common):
    def dependencies(self):
        return {'VE1': self._VE1}

    def calculate(self, VE1):
        return VE1 / self.mol.GetNumAtoms()


@method
class VE3(Common):
    def dependencies(self):
        return {'VE1': self._VE1}

    def calculate(self, VE1):
        if VE1 == 0:
            return MordredValueError('VE3: VE1 == zero')
        else:
            return np.log(0.1 * self.mol.GetNumAtoms() * VE1)


@method
class VR1(Common):
    def calculate(self, eig):
        s = 0.0

        for bond in self.mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            s += (eig.vec[i, eig.max] * eig.vec[j, eig.max]) ** -0.5

        return s


@method
class VR2(Common):
    def dependencies(self):
        return {'VR1': self._VR1}

    def calculate(self, VR1):
        return VR1 / self.mol.GetNumAtoms()


@method
class VR3(Common):
    def dependencies(self):
        return {'VR1': self._VR1}

    def calculate(self, VR1):
        if VR1 == 0:
            return MordredValueError('VR3: VR1 == 0')
        else:
            return np.log(0.1 * self.mol.GetNumAtoms() * VR1)


method_dict = {m.__name__: m for m in methods}


def get_method(n):
    if isinstance(n, string_types):
        return method_dict[n]

    return n
