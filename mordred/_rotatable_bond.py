import numpy as np

from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

from ._base import Descriptor
from ._bond_types import BondCount


class RotatableBondsBase(Descriptor):
    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()


class RotatableBondsCount(RotatableBondsBase):
    r"""ratatable bonds count descriptor(rdkit wrapper).

    :rtype: int
    """

    __slots__ = ()

    def __str__(self):
        return 'nRot'

    def calculate(self, mol):
        return CalcNumRotatableBonds(mol)


class RotatableBondsRatio(RotatableBondsBase):
    r"""rotatable bonds ratio descriptor.

    .. math::
        {\rm RotRatio} = \frac{N_{\rm rotatable bonds}}{N_{\rm bonds}}

    :rtype: float
    :returns: NaN when :math:`N_{\rm bonds} = 0`
    """

    __slots__ = ()

    def __str__(self):
        return 'RotRatio'

    def dependencies(self):
        return dict(
            nRot=RotatableBondsCount(),
            nB=BondCount('heavy'),
        )

    def calculate(self, mol, nRot, nB):
        if nB == 0:
            return np.nan

        return float(nRot) / float(nB)
