from ._base import Descriptor
from rdkit import Chem

class DistanceMatrix(Descriptor):
    @property
    def descriptor_key(self):
        return self.make_key(
            self.explicit_hydrogens,
            self.useBO, self.useAtomWts)

    def __init__(self, explicit_hydrogens, useBO, useAtomWts):
        self.explicit_hydrogens = explicit_hydrogens
        self.useBO = useBO
        self.useAtomWts = useAtomWts

    def calculate(self, mol):
        return Chem.GetDistanceMatrix(
            mol, useBO=self.useBO, useAtomWts=self.useAtomWts)
