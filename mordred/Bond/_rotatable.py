from .._base import Descriptor
from ._bond_types import BondCount

from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds


class RotatableBondsCount(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'nRot'

    def calculate(self, mol):
        return CalcNumRotatableBonds(mol)


class RotatableBondsRatio(Descriptor):
    explicit_hydrogens = False
    descriptor_name = 'RotRatio'

    @property
    def dependencies(self):
        return dict(
            nRot=RotatableBondsCount.make_key(),
            nB=BondCount.make_key('O'),
        )

    def calculate(self, mol, nRot, nB):
        return float(nRot) / float(nB)
