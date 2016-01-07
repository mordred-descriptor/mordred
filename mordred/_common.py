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


class AdjacencyMatrix(Descriptor):
    @property
    def descriptor_key(self):
        return self.make_key(
            self.explicit_hydrogens,
            self.useBO,
        )

    def __init__(self, explicit_hydrogens, useBO):
        self.explicit_hydrogens = explicit_hydrogens
        self.useBO = useBO

    def calculate(self, mol):
        return Chem.GetAdjacencyMatrix(mol, useBO=self.useBO)


class Valence(AdjacencyMatrix):
    @property
    def dependencies(self):
        return dict(D=AdjacencyMatrix.make_key(
            self.explicit_hydrogens,
            self.useBO,
        ))

    def calculate(self, mol, D):
        return D.sum(axis=0)
