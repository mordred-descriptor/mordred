from .._base import Descriptor
import numpy as np
from collections import namedtuple


Eig = namedtuple('eigen', ['val', 'vec'])

methods = []


def method(cls):
    methods.append(cls)
    return cls


class common(Descriptor):
    @property
    def explicitHydrogens(self):
        return self._explicitHydrogens

    @property
    def gasteigerCharges(self):
        return self._gasteigerCharges

    @property
    def kekulize(self):
        return self._kelulize

    @property
    def descriptor_key(self):
        return self.make_key(
            self.matrix,
            self.explicitHydrogens,
            self.gasteigerCharges,
            self.kekulize)

    def __init__(self, matrix, explicitHydrogens, gasteigerCharges, kekulize):
        self.matrix = matrix
        self._explicitHydrogens = explicitHydrogens
        self._gasteigerCharges = gasteigerCharges
        self._kelulize = kekulize

    @property
    def _key_args(self):
        return [self.matrix,
                self.explicitHydrogens,
                self.gasteigerCharges,
                self.kekulize]

    @property
    def _eig(self):
        return Eigen.make_key(*self._key_args)

    @property
    def _SpMax(self):
        return SpMax.make_key(*self._key_args)

    @property
    def _SpMean(self):
        return SpMean.make_key(*self._key_args)

    @property
    def _SpAD(self):
        return SpAD.make_key(*self._key_args)

    @property
    def _VE1(self):
        return VE1.make_key(*self._key_args)

    @property
    def _VR1(self):
        return VR1.make_key(*self._key_args)


class Eigen(common):
    @property
    def dependencies(self):
        return dict(matrix=self.matrix)

    def calculate(self, mol, matrix):
        w, v = np.linalg.eig(matrix)
        idx = w.argsort()
        w = w[idx]
        v = v[:, idx]
        if np.iscomplexobj(v):
            return Eig(w.real, v.real)
        else:
            return Eig(w, v)


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
        return eig.val[-1]


@method
class SpDiam(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig, SpMax=self._SpMax)

    def calculate(self, mol, SpMax, eig):
        return SpMax - eig.val[0]


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
class EE(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        return np.log(1+np.exp(eig.val).sum())


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
        return np.abs(eig.vec[:, 0].sum())


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
        try:
            return 0.1 * mol.GetNumAtoms() * np.log(VE1)
        except ValueError:
            return 0.0


@method
class VR1(common):
    @property
    def dependencies(self):
        return dict(eig=self._eig)

    def calculate(self, mol, eig):
        I, J = [], []
        for bond in mol.GetBonds():
            I.append(bond.GetBeginAtomIdx())
            J.append(bond.GetEndAtomIdx())

        m = np.abs(eig.vec[I, 0] * eig.vec[J, 0])
        return np.power(np.maximum(m, 1e-12), -0.5).sum()


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
        try:
            return 0.1 * mol.GetNumAtoms() * np.log(VR1)
        except ValueError:
            return 0.0
