from ._atomic_property import get_mc_gowan_volume
from ._base import Descriptor


class McGowanVolume(Descriptor):
    r"""McGowan volume descriptor.

    References
        * :cite:`10.1007/BF02311772`
    """

    __slots__ = ()

    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return 'VMcGowan'

    def __reduce_ex__(self, version):
        return self.__class__, ()

    def calculate(self, mol):
        a = sum(get_mc_gowan_volume(a) for a in mol.GetAtoms())
        return a - mol.GetNumBonds() * 6.56

    rtype = float
