from .._base import Descriptor
from ..Property import Weight, WildmanCrippenLogP

from rdkit.Chem import Lipinski as L


# Bioavailability: PSA, fusedAromaticRingCount
# LeadLikeness: LogD, ringCount
# MueggeFilter: ringCount, PSA
# VeberFilter: PSA


class LipinskiLike(Descriptor):
    '''
    Lipinski like descriptor

    Parameters:
        variant(str):
            * Lipinski
            * GooseFilter

    Returns:
        bool: filter result
    '''

    @property
    def descriptor_name(self):
        return self.variant

    @property
    def descriptor_key(self):
        return self.make_key(self.variant)

    def __init__(self, variant='Lipinski'):
        assert variant in set([
            'Lipinski',
            'Bioavailability',
            'GooseFilter',
            'LeadLikeness',
            'MueggeFilter',
            'VeberFilter',
        ])

        self.variant = variant

    @property
    def dependencies(self):
        return dict(
            LogP=WildmanCrippenLogP.make_key('LogP'),
            Weight=Weight.make_key(),
        )

    def Lipinski(self, mol, LogP, Weight, **other):
        return\
            L.NumHDonors(mol) <= 5 and\
            L.NumHAcceptors(mol) <= 10 and\
            Weight <= 500 and\
            LogP <= 5

    def calculate(self, mol, **deps):
        return getattr(self, self.variant)(mol, **deps)


_descriptors = [LipinskiLike]
__all__ = [d.__name__ for d in _descriptors]
