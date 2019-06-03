import numpy as np

from ._base import Descriptor
from ._util import atoms_to_numpy

__all__ = ("MomentOfInertia",)


class MomentOfInertiaBase(Descriptor):
    __slots__ = ()
    require_3D = True

    def _numpy(self):
        ws = atoms_to_numpy(lambda a: a.GetMass(), self.mol)
        ps = self.coord - np.sum(ws[:, np.newaxis] * self.coord, axis=0) / np.sum(ws)

        return ws, ps


class PrincipalAxis(MomentOfInertiaBase):
    __slots__ = ()

    def parameters(self):
        return ()

    def calculate(self):
        ws, ps = self._numpy()

        i = np.sum(
            -ws[:, np.newaxis, np.newaxis] * (ps[:, np.newaxis] * ps[:, :, np.newaxis]),
            axis=0,
        )

        diag = np.sum(
            ws[:, np.newaxis] * (np.sum(ps ** 2, axis=1)[:, np.newaxis] - ps ** 2),
            axis=0,
        )

        np.fill_diagonal(i, diag)
        return np.sort(np.linalg.eig(i)[0])[::-1]


class MomentOfInertia(MomentOfInertiaBase):
    since = "1.0.0"
    __slots__ = ("_axis",)

    def description(self):
        return "moment of inertia (axis = {})".format(self._axis)

    @classmethod
    def preset(cls, version):
        return map(cls, cls.axes)

    def __str__(self):
        return "MOMI-{}".format(self._axis)

    def parameters(self):
        return (self._axis,)

    axes = ("X", "Y", "Z")
    _axis_to_index = {a: i for i, a in enumerate(axes)}

    def __init__(self, axis="X"):
        assert axis in self.axes
        self._axis = axis

    def dependencies(self):
        return {"I": PrincipalAxis()}

    def calculate(self, I):
        return I[self._axis_to_index[self._axis]]

    rtype = float
