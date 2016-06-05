from ._base import Descriptor
from ._graph_matrix import Diameter3D as CDiameter3D
from ._graph_matrix import Radius3D as CRadius3D


__all__ = ('Diameter3D', 'Radius3D', 'GeometricalShapeIndex', 'PetitjeanIndex3D',)


class GeometricalIndexBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = True
    require_3D = True

    @classmethod
    def preset(cls):
        yield cls()

    def as_key(self):
        return self.__class__, ()

    rtype = float


class Radius3D(GeometricalIndexBase):
    r"""geometric radius descriptor."""
    __slots__ = ()

    def __str__(self):
        return 'GeomRadius'

    def dependencies(self):
        return {'R': CRadius3D(self.explicit_hydrogens)}

    def calculate(self, R):
        return R


class Diameter3D(GeometricalIndexBase):
    r"""geometric diameter descriptor."""
    __slots__ = ()

    def __str__(self):
        return 'GeomDiameter'

    def dependencies(self):
        return {'D': CDiameter3D(self.explicit_hydrogens)}

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
    __slots__ = ()

    def __str__(self):
        return 'GeomShapeIndex'

    def dependencies(self):
        return {
            'R': CRadius3D(self.explicit_hydrogens),
            'D': CDiameter3D(self.explicit_hydrogens),
        }

    def calculate(self, R, D):
        with self.rethrow_zerodiv():
            return float(D - R) / float(R)


class PetitjeanIndex3D(GeometricalShapeIndex):
    r"""geometric Petitjean index descriptor.

    .. math::

        I_{\rm Petitjean} = \frac{D - R}{D}

    where
    :math:`R` is geometric radius,
    :math:`D` is geometric diameter.

    :returns: NaN when :math:`D = 0`
    """
    __slots__ = ()

    def __str__(self):
        return 'GeomPetitjeanIndex'

    def calculate(self, R, D):
        with self.rethrow_zerodiv():
            return float(D - R) / float(D)
