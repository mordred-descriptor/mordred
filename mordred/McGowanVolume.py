from ._base import Descriptor
from ._atomic_property import get_mc_gowan_volume

__all__ = ("McGowanVolume",)


class McGowanVolume(Descriptor):
    r"""McGowan volume descriptor.

    References
        * :doi:`10.1007/BF02311772`

    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "McGowan volume"

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return "VMcGowan"

    def parameters(self):
        return ()

    def calculate(self):
        a = sum(get_mc_gowan_volume(a) for a in self.mol.GetAtoms())
        return a - self.mol.GetNumBonds() * 6.56

    rtype = float
