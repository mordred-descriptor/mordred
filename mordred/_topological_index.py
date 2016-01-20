from ._base import Descriptor
from ._common import Diameter as CDiameter
from ._common import Radius as CRadius


class TopologicalIndexBase(Descriptor):
    explicit_hydrogens = False
    require_connected = False

    @classmethod
    def preset(cls):
        yield cls()


class Radius(TopologicalIndexBase):
    r"""radius descriptor.

    :rtype: int
    """

    def __str__(self):
        return 'Radius'

    def dependencies(self):
        return dict(
            R=CRadius(
                self.explicit_hydrogens,
                False, False,
            )
        )

    def calculate(self, mol, R):
        return int(R)


class Diameter(TopologicalIndexBase):
    r"""diameter descriptor.

    :rtype: int
    """

    def __str__(self):
        return 'Diameter'

    def dependencies(self):
        return dict(
            D=CDiameter(
                self.explicit_hydrogens,
                False, False)
        )

    def calculate(self, mol, D):
        return int(D)


class TopologicalShapeIndex(TopologicalIndexBase):
    r"""topological shape index descriptor.

    .. math::

        I_{\rm topo} = \frac{D - R}{R}

    where
    :math:`R` is graph radius,
    :math:`D` is graph diameter.

    :rtype: float
    """

    def __str__(self):
        return 'TopoShapeIndex'

    def dependencies(self):
        args = [self.explicit_hydrogens, False, False]

        return dict(
            R=CRadius(*args),
            D=CDiameter(*args)
        )

    def calculate(self, mol, R, D):
        return float(D - R) / float(R)


class PetitjeanIndex(TopologicalIndexBase):
    r"""Petitjean index descriptor.

    .. math::

        I_{\rm Petitjean} = \frac{D - R}{D}

    where
    :math:`R` is graph radius,
    :math:`D` is graph diameter.

    :rtype: float
    """

    def __str__(self):
        return 'PetitjeanIndex'

    def dependencies(self):
        args = [self.explicit_hydrogens, False, False]

        return dict(
            R=CRadius(*args),
            D=CDiameter(*args)
        )

    def calculate(self, mol, R, D):
        return float(D - R) / float(D)
