import numpy as np

from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

from .BondCount import BondCount

from ._base import Descriptor


__all__ = (
    'RotatableBondsBase', 'RotatableBondsCount', 'RotatableBondsRatio',
)


class RotatableBondsBase(Descriptor):
    explicit_hydrogens = False

    @classmethod
    def preset(cls):
        yield cls()

    def as_key(self):
        return self.__class__, ()


class RotatableBondsCount(RotatableBondsBase):
    r"""ratatable bonds count descriptor(rdkit wrapper)."""

    __slots__ = ()

    def __str__(self):
        return 'nRot'

    def calculate(self):
        return CalcNumRotatableBonds(self.mol)

    rtype = int


class RotatableBondsRatio(RotatableBondsBase):
    r"""rotatable bonds ratio descriptor.

    .. math::
        {\rm RotRatio} = \frac{N_{\rm rotatable bonds}}{N_{\rm bonds}}

    :returns: NaN when :math:`N_{\rm bonds} = 0`
    """

    __slots__ = ()

    def __str__(self):
        return 'RotRatio'

    def dependencies(self):
        return {
            'nRot': RotatableBondsCount(),
            'nB': BondCount('heavy'),
        }

    def calculate(self, nRot, nB):
        if nB == 0:
            return np.nan

        return float(nRot) / float(nB)

    rtype = float
