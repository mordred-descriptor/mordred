import numpy as np

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
            mol, useBO=self.useBO, useAtomWts=self.useAtomWts
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


class Diameter(DistanceMatrix):
    def dependencies(self):
        return dict(
            D=DistanceMatrix(
                self.explicit_hydrogens,
                self.useBO,
                self.useAtomWts,
            )
        )

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
            return Chem.GetAdjacencyMatrix(mol, useBO=self.useBO)

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
