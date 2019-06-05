from rdkit.Chem import MolSurf
from rdkit.Chem.EState import EState_VSA as RDKit_EState_VSA

from ._base import Descriptor

__all__ = ("LabuteASA", "PEOE_VSA", "SMR_VSA", "SlogP_VSA", "EState_VSA", "VSA_EState")


class LabuteASA(Descriptor):
    r"""Labute's Approximate Surface Area descriptor(rdkit wrapper)."""

    since = "1.0.0"
    __slots__ = ()
    explicit_hydrogens = False

    def description(self):
        return "Labute's Approximate Surface Area"

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return self.__class__.__name__

    def parameters(self):
        return ()

    def calculate(self):
        return MolSurf.LabuteASA(self.mol)

    rtype = float


class MoeTypeBase(Descriptor):
    __slots__ = ("_k",)
    explicit_hydrogens = False
    _module = MolSurf

    @classmethod
    def preset(cls, version):
        return map(cls, range(1, cls.k_max))

    def description(self):
        return self._fn.__doc__

    @property
    def _fn(self):
        return getattr(self._module, str(self))

    def __str__(self):
        return self.__class__.__name__ + str(self._k)

    def parameters(self):
        return (self._k,)

    def __init__(self, k=1):
        assert 1 <= k <= self.k_max
        self._k = k

    def calculate(self):
        return self._fn(self.mol)

    rtype = float


class PEOE_VSA(MoeTypeBase):
    r"""MOE type descriptors using gasteiger charge and surface area contribution(rdkit wrapper).

    :type k: int
    :param k: (:math:`1 <= k <= k_{\rm max}`)
    """

    since = "1.0.0"
    __slots__ = ()
    k_max = 14


class SMR_VSA(MoeTypeBase):
    r"""MOE type descriptors using Wildman-Crippen MR and surface area contribution(rdkit wrapper).

    :type k: int
    :param k: (:math:`1 <= k <= k_{\rm max}`)
    """

    since = "1.0.0"
    __slots__ = ()
    k_max = 10


class SlogP_VSA(MoeTypeBase):
    r"""MOE type descriptors using Wildman-Crippen LogP and surface area contribution(rdkit wrapper).

    :type k: int
    :param k: (:math:`1 <= k <= k_{\rm max}`)
    """

    since = "1.0.0"
    __slots__ = ()
    k_max = 12


class EState_VSA(MoeTypeBase):
    r"""MOE type descriptors using EState indices and surface area contribution(rdkit wrapper).

    :type k: int
    :param k: (:math:`1 <= k <= k_{\rm max}`)
    """

    since = "1.0.0"
    __slots__ = ()
    _module = RDKit_EState_VSA
    k_max = 11


class VSA_EState(MoeTypeBase):
    r"""MOE type descriptors using EState indices and surface area contribution(rdkit wrapper).

    :type k: int
    :param k: (:math:`1 <= k <= k_{\rm max}`)
    """

    since = "1.0.0"
    __slots__ = ()
    _module = RDKit_EState_VSA
    k_max = 10
