from ._base import Descriptor
from .Bond import HBondDonor, HBondAcceptor
from .Property import Weight, WildmanCrippenLogP


class LipinskiLike(Descriptor):
    r'''
    Lipinski like descriptor

    LogP: WildmanCrippenLogP

    :type variant: str
    :param variant:
        * 'Lipinski'
        * 'GhoseFilter'

    :rtype: bool
    '''

    require_connected = False

    @classmethod
    def preset(cls):
        yield cls('Lipinski')
        yield cls('GhoseFilter')

    def __str__(self):
        return self.variant

    descriptor_keys = 'variant',

    def __init__(self, variant='Lipinski'):
        assert variant in set([
            'Lipinski',
            'GhoseFilter',
        ])

        self.variant = variant

    def dependencies(self):
        return {
            prop: key
            for prop, key in deps_keys.items()
            if prop in deps_dict[self.variant]
        }

    def _Lipinski(self, mol, LogP, MW, HBDon, HBAcc):
        return\
            HBDon <= 5 and\
            HBAcc <= 10 and\
            MW <= 500 and\
            LogP <= 5

    def _GhoseFilter(self, mol, MW, LogP, MR):
        return\
            (160 <= MW <= 480) and\
            (20 <= mol.GetNumAtoms() <= 70) and\
            (-0.4 <= LogP <= 5.6) and\
            (40 <= MR <= 130)

    def calculate(self, mol, **deps):
        return getattr(self, '_' + self.variant)(mol, **deps)

deps_dict = dict(
    Lipinski=set(['LogP', 'MW', 'HBDon', 'HBAcc']),
    GhoseFilter=set(['MW', 'LogP', 'MR']),
)

deps_keys = dict(
    LogP=WildmanCrippenLogP('LogP'),
    MR=WildmanCrippenLogP('MR'),

    MW=Weight(),

    HBDon=HBondDonor(),
    HBAcc=HBondAcceptor(),
)
