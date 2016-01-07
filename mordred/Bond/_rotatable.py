from .._base import Descriptor
from ._bond_types import BondCount

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


class RotatableBondsRatio(Descriptor):
    explicit_hydrogens = False

    @property
    def descriptor_name(self):
        return 'RotRatio'

    @property
    def descriptor_key(self):
        return self.make_key()

    @property
    def dependencies(self):
        return dict(
            nRot=RotatableBondsCount.make_key(),
            nB=BondCount.make_key('O'),
        )

    def calculate(self, mol, nRot, nB):
        return float(nRot) / float(nB)
