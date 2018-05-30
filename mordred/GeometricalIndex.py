from __future__ import division

from ._base import Descriptor
from ._graph_matrix import Radius3D as CRadius3D
from ._graph_matrix import Diameter3D as CDiameter3D

__all__ = ("Diameter3D", "Radius3D", "GeometricalShapeIndex", "PetitjeanIndex3D")


class GeometricalIndexBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = True
    require_3D = True

    @classmethod
    def preset(cls, version):
        yield cls()

    def parameters(self):
        return ()

    rtype = float


class Radius3D(GeometricalIndexBase):
    r"""geometric radius descriptor."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "geometric radius"

    def __str__(self):
        return "GeomRadius"

    def dependencies(self):
        return {"R": CRadius3D(self.explicit_hydrogens)}

    def calculate(self, R):
        return R


class Diameter3D(GeometricalIndexBase):
    r"""geometric diameter descriptor."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "geometric diameter"

    def __str__(self):
        return "GeomDiameter"

    def dependencies(self):
        return {"D": CDiameter3D(self.explicit_hydrogens)}

    def calculate(self, D):
        return D


class GeometricalShapeIndex(GeometricalIndexBase):
    r"""geometrical shape index descriptor.

    .. math::

        I_{\rm topo} = \frac{D - R}{R}

    where
    :math:`R` is geometric radius,
    :math:`D` is geometric diameter.

    :returns: NaN when :math:`R = 0`
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "geometrical shape index"

    def __str__(self):
        return "GeomShapeIndex"

    def dependencies(self):
        return {
            "D": CDiameter3D(self.explicit_hydrogens),
            "R": CRadius3D(self.explicit_hydrogens),
        }

    def calculate(self, R, D):
        with self.rethrow_zerodiv():
            return (D - R) / (R)


class PetitjeanIndex3D(GeometricalShapeIndex):
    r"""geometric Petitjean index descriptor.

    .. math::

        I_{\rm Petitjean} = \frac{D - R}{D}

    where
    :math:`R` is geometric radius,
    :math:`D` is geometric diameter.

    :returns: NaN when :math:`D = 0`
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "geometric Petitjean index"

    def __str__(self):
        return "GeomPetitjeanIndex"

    def calculate(self, R, D):
        with self.rethrow_zerodiv():
            return (D - R) / (D)
