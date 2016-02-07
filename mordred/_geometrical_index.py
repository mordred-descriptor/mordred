from ._base import Descriptor
from ._common import Diameter3D as CDiameter3D
from ._common import Radius3D as CRadius3D


class GeometricalIndexBase(Descriptor):
    explicit_hydrogens = True
    require_3D = True

    @classmethod
    def preset(cls):
        yield cls()

    def __reduce_ex__(self, version):
        return self.__class__, ()

    rtype = float


class Radius3D(GeometricalIndexBase):
    r"""geometric radius descriptor."""

    __slots__ = ()

    def __str__(self):
        return 'GeomRadius'

    def dependencies(self):
        return dict(
            R=CRadius3D(self.explicit_hydrogens)
        )

    def calculate(self, mol, conf, R):
        return R


class Diameter3D(GeometricalIndexBase):
    r"""geometric diameter descriptor."""

    __slots__ = ()

    def __str__(self):
        return 'GeomDiameter'

    def dependencies(self):
        return dict(
            D=CDiameter3D(self.explicit_hydrogens)
        )

    def calculate(self, mol, conf, D):
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
        return dict(
            R=CRadius3D(self.explicit_hydrogens),
            D=CDiameter3D(self.explicit_hydrogens),
        )

    def calculate(self, mol, conf, R, D):
        if R == 0:
            return float('nan')

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

    def calculate(self, mol, conf, R, D):
        if D == 0:
            return float('nan')

        return float(D - R) / float(D)
