import numpy as np

from ._base import Descriptor
from ._util import atoms_to_numpy


__all__ = ('MomentOfInertia',)


class MomentOfInertiaBase(Descriptor):
    require_3D = True

    @staticmethod
    def _numpy(mol, conf):
        ws = atoms_to_numpy(lambda a: a.GetMass(), mol)
        ps = conf - np.sum(ws[:, np.newaxis] * conf, axis=0) / np.sum(ws)

        return ws, ps


class PrincipalAxis(MomentOfInertiaBase):
    def __reduce_ex__(self, version):
        return self.__class__, ()

    def calculate(self, mol, conf):
        ws, ps = self._numpy(mol, conf)

        I = np.sum(
            -ws[:, np.newaxis, np.newaxis] *
            (ps[:, np.newaxis] * ps[:, :, np.newaxis]),
            axis=0
        )

        diag = np.sum(
            ws[:, np.newaxis] *
            (np.sum(ps ** 2, axis=1)[:, np.newaxis] - ps ** 2),
            axis=0
        )

        np.fill_diagonal(I, diag)
        return np.sort(np.linalg.eig(I)[0])[::-1]


class MomentOfInertia(MomentOfInertiaBase):
    @classmethod
    def preset(cls):
        return map(cls, cls.axes)

    def __str__(self):
        return 'MOMI-{}'.format(self._axis)

    def __reduce_ex__(self, version):
        return self.__class__, (self._axis,)

    axes = ('X', 'Y', 'Z')
    _axis_to_index = {a: i for i, a in enumerate(axes)}

    def __init__(self, axis='X'):
        assert axis in self.axes
        self._axis = axis

    def dependencies(self):
        return {'I': PrincipalAxis()}

    def calculate(self, mol, conf, I):
        return I[self._axis_to_index[self._axis]]

    rtype = float


if __name__ == '__main__':
    from .__main__ import submodule
    submodule()
