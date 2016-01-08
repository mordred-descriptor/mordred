from .._base import Descriptor
from .._common import Radius as _R, Diameter as _D


class Radius(Descriptor):
    explicit_hydrogens = False

    @property
    def dependencies(self):
        return dict(
            R=_R.make_key(
                self.explicit_hydrogens,
                False, False)
        )

    @property
    def descriptor_name(self):
        return 'Radius'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol, R):
        return R


class Diameter(Descriptor):
    explicit_hydrogens = False

    @property
    def dependencies(self):
        return dict(
            D=_D.make_key(
                self.explicit_hydrogens,
                False, False)
        )

    @property
    def descriptor_name(self):
        return 'Diameter'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol, D):
        return D


class TopologicalShapeIndex(Descriptor):
    explicit_hydrogens = False

    @property
    def dependencies(self):
        args = [self.explicit_hydrogens, False, False]

        return dict(
            R=_R.make_key(*args),
            D=_D.make_key(*args)
        )

    @property
    def descriptor_name(self):
        return 'TopoShapeIndex'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol, R, D):
        return float(D - R) / float(R)
