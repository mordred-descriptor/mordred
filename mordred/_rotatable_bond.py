from ._base import Descriptor
from ._bond_types import BondCount

from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds


class RotatableBondsBase(Descriptor):
    explicit_hydrogens = False
    require_connected = False


class RotatableBondsCount(RotatableBondsBase):
    r'''
    ratatable bonds count descriptor

    Returns:
        int: rotatable bonds count
    '''

    def __str__(self):
        return 'nRot'

    def calculate(self, mol):
        return CalcNumRotatableBonds(mol)


class RotatableBondsRatio(RotatableBondsBase):
    r'''
    rotatable bonds ratio descriptor

    .. math::
        {\rm RotRatio} = \frac{N_{\rm rotatable bonds}}{N_{\rm bonds}}

    Returns:
        float: rotatable bonds ratio
    '''

    def __str__(self):
        return 'RotRatio'

    @property
    def dependencies(self):
        return dict(
            nRot=RotatableBondsCount(),
            nB=BondCount('heavy'),
        )

    def calculate(self, mol, nRot, nB):
        return float(nRot) / float(nB)
