from .._base import Descriptor
from .._atomic_property import get_mc_gowan_volume


class McGowanVolume(Descriptor):
    descriptor_name = 'McGowan'

    def calculate(self, mol):
        a = sum(get_mc_gowan_volume(a) for a in mol.GetAtoms())
        return a - mol.GetNumBonds() * 6.56
