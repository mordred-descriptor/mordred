from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

from ._base import Descriptor
from .BondCount import BondCount

__all__ = ("RotatableBondsCount", "RotatableBondsRatio")


class RotatableBondsBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False

    @classmethod
    def preset(cls, version):
        yield cls()

    def parameters(self):
        return ()


class RotatableBondsCount(RotatableBondsBase):
    r"""rotatable bonds count descriptor(rdkit wrapper)."""

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "rotatable bonds count"

    def __str__(self):
        return "nRot"

    def calculate(self):
        return CalcNumRotatableBonds(self.mol)

    rtype = int


class RotatableBondsRatio(RotatableBondsBase):
    r"""rotatable bonds ratio descriptor.

    .. math::
        {\rm RotRatio} = \frac{N_{\rm rotatable bonds}}{N_{\rm bonds}}

    :returns: NaN when :math:`N_{\rm bonds} = 0`
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "rotatable bonds ratio"

    def __str__(self):
        return "RotRatio"

    def dependencies(self):
        return {"nB": BondCount("heavy"), "nRot": RotatableBondsCount()}

    def calculate(self, nRot, nB):
        with self.rethrow_zerodiv():
            return float(nRot) / float(nB)

    rtype = float
