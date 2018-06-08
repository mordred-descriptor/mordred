from __future__ import division

from ._base import Descriptor
from .RingCount import Rings
from .molecular_framework import _get_linkers

__all__ = ("Framework",)


class FrameworkCache(Descriptor):
    __slots__ = ()

    def parameters(self):
        return ()

    def dependencies(self):
        return {"Rs": Rings()}

    def calculate(self, Rs):
        return _get_linkers(self.mol, Rs), Rs


class Framework(Descriptor):
    r"""molecular framework ratio descriptor.

    .. math::

        f_{\rm MF} = \frac{N_{\rm MF}}{N}

    where
    :math:`N_{\rm MF}` is number of atoms in molecular framework,
    :math:`N` is number of all atoms.

    References
        * :doi:`10.1021/jm9602928`

    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "molecular framework ratio"

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return "fMF"

    def parameters(self):
        return ()

    def dependencies(self):
        return {"F": FrameworkCache()}

    def calculate(self, F):
        linkers, rings = F
        indices = {i for linker in linkers for ab in linker for i in ab}
        indices.update(i for ring in rings for i in ring)
        N = self.mol.GetNumAtoms()

        return len(indices) / N

    rtype = float
