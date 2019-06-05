import numpy as np
from rdkit import Chem

from ._base import Descriptor


class DistanceMatrix(Descriptor):
    __slots__ = ("explicit_hydrogens", "useBO", "useAtomWts")

    hermitian = True

    def parameters(self):
        return self.explicit_hydrogens, self.useBO, self.useAtomWts

    def __init__(self, explicit_hydrogens, useBO=False, useAtomWts=False):
        self.explicit_hydrogens = explicit_hydrogens
        self.useBO = useBO
        self.useAtomWts = useAtomWts

    def calculate(self):
        return Chem.GetDistanceMatrix(
            self.mol, useBO=self.useBO, useAtomWts=self.useAtomWts, force=True
        )


class Eccentricity(DistanceMatrix):
    __slots__ = ()

    def dependencies(self):
        return {
            "D": DistanceMatrix(self.explicit_hydrogens, self.useBO, self.useAtomWts)
        }

    def calculate(self, D):
        return D.max(axis=0)


class Radius(Eccentricity):
    __slots__ = ()

    def dependencies(self):
        return {"E": Eccentricity(self.explicit_hydrogens, self.useBO, self.useAtomWts)}

    def calculate(self, E):
        return E.min()


class Diameter(Eccentricity):
    __slots__ = ()

    def calculate(self, D):
        return D.max()


class AdjacencyMatrix(Descriptor):
    __slots__ = ("explicit_hydrogens", "useBO", "order")

    hermitian = True

    def parameters(self):
        return self.explicit_hydrogens, self.useBO, self.order

    def __init__(self, explicit_hydrogens, useBO=False, order=1):
        self.explicit_hydrogens = explicit_hydrogens
        self.useBO = useBO
        self.order = order

    def dependencies(self):
        if self.order > 1:
            return {
                "A1": self.__class__(self.explicit_hydrogens, self.useBO, 1),
                "An": self.__class__(
                    self.explicit_hydrogens, self.useBO, self.order - 1
                ),
            }
        else:
            return {}

    def calculate(self, An=None, A1=None):
        if self.order == 1:
            return Chem.GetAdjacencyMatrix(self.mol, useBO=self.useBO, force=True)

        return An.dot(A1)


class Valence(AdjacencyMatrix):
    __slots__ = ()

    def dependencies(self):
        return {"D": AdjacencyMatrix(self.explicit_hydrogens, self.useBO)}

    def calculate(self, D):
        return D.sum(axis=0)


class DistanceMatrix3D(Descriptor):
    __slots__ = ("explicit_hydrogens", "useAtomWts")

    require_3D = True

    def parameters(self):
        return self.explicit_hydrogens, self.useAtomWts

    def __init__(self, explicit_hydrogens, useAtomWts=False):
        self.explicit_hydrogens = explicit_hydrogens
        self.useAtomWts = useAtomWts

    def calculate(self):
        return np.sqrt(np.sum((self.coord[:, np.newaxis] - self.coord) ** 2, axis=2))


class Eccentricity3D(DistanceMatrix3D):
    __slots__ = ()

    def dependencies(self):
        return {"D": DistanceMatrix3D(self.explicit_hydrogens, self.useAtomWts)}

    def calculate(self, D):
        return D.max(axis=0)


class Radius3D(Eccentricity3D):
    __slots__ = ()

    def dependencies(self):
        return {"E": Eccentricity3D(self.explicit_hydrogens, self.useAtomWts)}

    def calculate(self, E):
        return E.min()


class Diameter3D(Eccentricity3D):
    __slots__ = ()

    def calculate(self, D):
        return D.max()
