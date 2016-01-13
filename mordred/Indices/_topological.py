from .._base import Descriptor
from .._common import Radius as _R, Diameter as _D


class Radius(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'Radius'

    @property
    def dependencies(self):
        return dict(
            R=_R.make_key(
                self.explicit_hydrogens,
                False, False)
        )

    def calculate(self, mol, R):
        return int(R)


class Diameter(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'Diameter'

    @property
    def dependencies(self):
        return dict(
            D=_D.make_key(
                self.explicit_hydrogens,
                False, False)
        )

    def calculate(self, mol, D):
        return int(D)


class TopologicalShapeIndex(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'TopoShapeIndex'

    @property
    def dependencies(self):
        args = [self.explicit_hydrogens, False, False]

        return dict(
            R=_R.make_key(*args),
            D=_D.make_key(*args)
        )

    def calculate(self, mol, R, D):
        return float(D - R) / float(R)


class PetitjeanIndex(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'PetitjeanIndex'

    @property
    def dependencies(self):
        args = [self.explicit_hydrogens, False, False]

        return dict(
            R=_R.make_key(*args),
            D=_D.make_key(*args)
        )

    def calculate(self, mol, R, D):
        return float(D - R) / float(D)
