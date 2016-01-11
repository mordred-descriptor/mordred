from ._base import Descriptor
from rdkit import Chem


class DistanceMatrix(Descriptor):
    @property
    def descriptor_key(self):
        return self.make_key(
            self.explicit_hydrogens,
            self.useBO,
            self.useAtomWts,
        )

    def __init__(self, explicit_hydrogens, useBO, useAtomWts):
        self.explicit_hydrogens = explicit_hydrogens
        self.useBO = useBO
        self.useAtomWts = useAtomWts

    def calculate(self, mol):
        return Chem.GetDistanceMatrix(
            mol, useBO=self.useBO, useAtomWts=self.useAtomWts)


class Eccentricity(DistanceMatrix):
    @property
    def dependencies(self):
        return dict(D=DistanceMatrix.make_key(
            self.explicit_hydrogens,
            self.useBO,
            self.useAtomWts,
        ))

    def calculate(self, mol, D):
        return D.max(axis=0)


class Radius(Eccentricity):
    @property
    def dependencies(self):
        return dict(E=Eccentricity.make_key(
            self.explicit_hydrogens,
            self.useBO,
            self.useAtomWts,
        ))

    def calculate(self, mol, E):
        return E.min()


class Diameter(DistanceMatrix):
    @property
    def dependencies(self):
        return dict(D=DistanceMatrix.make_key(
            self.explicit_hydrogens,
            self.useBO,
            self.useAtomWts,
        ))

    def calculate(self, mol, D):
        return D.max()


class AdjacencyMatrix(Descriptor):
    @property
    def descriptor_key(self):
        return self.make_key(
            self.explicit_hydrogens,
            self.useBO,
            self.order,
        )

    def __init__(self, explicit_hydrogens, useBO, order=1):
        self.explicit_hydrogens = explicit_hydrogens
        self.useBO = useBO
        self.order = order

    @property
    def dependencies(self):
        if self.order > 1:
            return dict(
                An=self.make_key(
                    self.explicit_hydrogens,
                    self.useBO,
                    self.order - 1,
                ),
                A1=self.make_key(
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
    @property
    def dependencies(self):
        return dict(D=AdjacencyMatrix.make_key(
            self.explicit_hydrogens,
            self.useBO,
        ))

    def calculate(self, mol, D):
        return D.sum(axis=0)
