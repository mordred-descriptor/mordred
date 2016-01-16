from ._base import Descriptor
from ._atomic_property import get_mc_gowan_volume


class McGowanVolume(Descriptor):
    r'''
    McGowan volume descriptor

    Returns:
        float: McGowan volume
    '''

    def __str__(self):
        return 'VMcGowan'

    def calculate(self, mol):
        a = sum(get_mc_gowan_volume(a) for a in mol.GetAtoms())
        return a - mol.GetNumBonds() * 6.56
