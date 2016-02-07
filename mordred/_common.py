from rdkit import Chem

from ._base import Descriptor


class DistanceMatrix(Descriptor):
    __slots__ = ('explicit_hydrogens', 'useBO', 'useAtomWts',)

    def __reduce_ex__(self, version):
        return self.__class__, (
            self.explicit_hydrogens,
            self.useBO,
            self.useAtomWts,
        )

    def __init__(self, explicit_hydrogens, useBO=False, useAtomWts=False):
        self.explicit_hydrogens = explicit_hydrogens
        self.useBO = useBO
        self.useAtomWts = useAtomWts

    def calculate(self, mol):
        return Chem.GetDistanceMatrix(
            mol, useBO=self.useBO, useAtomWts=self.useAtomWts, force=True,
        )


class Eccentricity(DistanceMatrix):
    def dependencies(self):
        return dict(
            D=DistanceMatrix(
                self.explicit_hydrogens,
                self.useBO,
                self.useAtomWts,
            )
        )

    def calculate(self, mol, D):
        return D.max(axis=0)


class Radius(Eccentricity):
    def dependencies(self):
        return dict(
            E=Eccentricity(
                self.explicit_hydrogens,
                self.useBO,
                self.useAtomWts,
            )
        )

    def calculate(self, mol, E):
        return E.min()


class Diameter(Eccentricity):
    def calculate(self, mol, D):
        return D.max()


class AdjacencyMatrix(Descriptor):
    __slots__ = ('explicit_hydrogens', 'useBO', 'order',)

    def __reduce_ex__(self, version):
        return self.__class__, (
            self.explicit_hydrogens,
            self.useBO,
            self.order,
        )

    def __init__(self, explicit_hydrogens, useBO=False, order=1):
        self.explicit_hydrogens = explicit_hydrogens
        self.useBO = useBO
        self.order = order

    def dependencies(self):
        if self.order > 1:
            return dict(
                An=self.__class__(
                    self.explicit_hydrogens,
                    self.useBO,
                    self.order - 1,
                ),
                A1=self.__class__(
                    self.explicit_hydrogens,
                    self.useBO,
                    1
                )
            )
        else:
            return dict()

    def calculate(self, mol, An=None, A1=None):
        if self.order == 1:
            return Chem.GetAdjacencyMatrix(mol, useBO=self.useBO, force=True)

        return An.dot(A1)


class Valence(AdjacencyMatrix):
    def dependencies(self):
        return dict(
            D=AdjacencyMatrix(
                self.explicit_hydrogens,
                self.useBO,
            )
        )

    def calculate(self, mol, D):
        return D.sum(axis=0)


class DistanceMatrix3D(Descriptor):
    __slots__ = ('explicit_hydrogens', 'useAtomWts',)

    require_3D = True

    def __reduce_ex__(self, version):
        return self.__class__, (
            self.explicit_hydrogens,
            self.useAtomWts,
        )

    def __init__(self, explicit_hydrogens, useAtomWts=False):
        self.explicit_hydrogens = explicit_hydrogens
        self.useAtomWts = useAtomWts

    def calculate(self, mol, conf):
        return Chem.Get3DDistanceMatrix(
            mol, confId=conf.GetId(),
            useAtomWts=self.useAtomWts, force=True,
        )


class Eccentricity3D(DistanceMatrix3D):
    def dependencies(self):
        return dict(
            D=DistanceMatrix3D(
                self.explicit_hydrogens,
                self.useAtomWts,
            )
        )

    def calculate(self, mol, conf, D):
        return D.max(axis=0)


class Radius3D(Eccentricity3D):
    def dependencies(self):
        return dict(
            E=Eccentricity3D(
                self.explicit_hydrogens,
                self.useAtomWts,
            )
        )

    def calculate(self, mol, conf, E):
        return E.min()


class Diameter3D(Eccentricity3D):
    def calculate(self, mol, conf, D):
        return D.max()
