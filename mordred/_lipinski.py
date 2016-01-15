from ._base import Descriptor
from .Bond import HBondDonor, HBondAcceptor
from .Property import Weight, WildmanCrippenLogP
from .TPSA import TPSA

from rdkit.Chem import Lipinski as L


# Bioavailability: PSA, fusedAromaticRingCount
# LeadLikeness: LogD, ringCount
# MueggeFilter: ringCount, PSA
# VeberFilter: PSA


class LipinskiLike(Descriptor):
    '''
    Lipinski like descriptor

    LogP: WildmanCrippenLogP

    Parameters:
        variant(str):
            * Lipinski
            * GhoseFilter

    Returns:
        bool: filter result
    '''

    @property
    def descriptor_name(self):
        return self.variant

    @property
    def descriptor_key(self):
        return self.make_key(self.variant)

    @classmethod
    def preset(cls):
        yield cls('Lipinski')
        yield cls('GhoseFilter')

    def __init__(self, variant='Lipinski'):
        assert variant in set([
            'Lipinski',
            'GhoseFilter',
        ])

        self.variant = variant

    @property
    def dependencies(self):
        return {
            prop: key
            for prop, key in deps_keys.items()
            if prop in deps_dict[self.variant]
        }


    def Lipinski(self, mol, LogP, MW, HBDon, HBAcc):
        return\
            HBDon <= 5 and\
            HBAcc <= 10 and\
            MW <= 500 and\
            LogP <= 5

    def GhoseFilter(self, mol, MW, LogP, MR):
        return\
            (160 <= MW <= 480) and\
            (20 <= mol.GetNumAtoms() <= 70) and\
            (-0.4 <= LogP <= 5.6) and\
            (40 <= MR <= 130)

    def calculate(self, mol, **deps):
        return getattr(self, self.variant)(mol, **deps)

deps_dict = dict(
    Lipinski=set(['LogP', 'MW', 'HBDon', 'HBAcc']),
    GhoseFilter=set(['MW', 'LogP', 'MR']),
)

deps_keys = dict(
    LogP=WildmanCrippenLogP.make_key('LogP'),
    MR=WildmanCrippenLogP.make_key('MR'),

    MW=Weight.make_key(),

    HBDon=HBondDonor.make_key(),
    HBAcc=HBondAcceptor.make_key(),
)
