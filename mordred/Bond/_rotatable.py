from .._base import Descriptor
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds


class RotatableBondsCount(Descriptor):
    explicit_hydrogens = False

    @property
    def descriptor_name(self):
        return 'nRot'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol):
        return CalcNumRotatableBonds(mol)
