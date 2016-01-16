from ._base import Descriptor
from rdkit.Chem import EState as ES
from numpy import nan
from functools import reduce

try:
    import builtins
except ImportError:
    import __builtin__ as builtins


es_types = [
    'sLi', 'ssBe', 'ssssBe', 'ssBH', 'sssB', 'ssssB', 'sCH3', 'dCH2', 'ssCH2',
    'tCH', 'dsCH', 'aaCH', 'sssCH', 'ddC', 'tsC', 'dssC', 'aasC', 'aaaC',
    'ssssC', 'sNH3', 'sNH2', 'ssNH2', 'dNH', 'ssNH', 'aaNH', 'tN', 'sssNH',
    'dsN', 'aaN', 'sssN', 'ddsN', 'aasN', 'ssssN', 'sOH', 'dO', 'ssO', 'aaO',
    'sF', 'sSiH3', 'ssSiH2', 'sssSiH', 'ssssSi', 'sPH2', 'ssPH', 'sssP',
    'dsssP', 'sssssP', 'sSH', 'dS', 'ssS', 'aaS', 'dssS', 'ddssS', 'sCl',
    'sGeH3', 'ssGeH2', 'sssGeH', 'ssssGe', 'sAsH2', 'ssAsH', 'sssAs', 'sssdAs',
    'sssssAs', 'sSeH', 'dSe', 'ssSe', 'aaSe', 'dssSe', 'ddssSe', 'sBr', 'sSnH3',
    'ssSnH2', 'sssSnH', 'ssssSn', 'sI', 'sPbH3', 'ssPbH2', 'sssPbH', 'ssssPb'
]

es_type_set = set(es_types)


class EStateBase(Descriptor):
    explicit_hydrogens = False


class EStateCache(EStateBase):
    def calculate(self, mol):
        return ES.TypeAtoms(mol), ES.EStateIndices(mol)


aggr_name_dict = dict(
    count='N',
    sum='S',
    max='MAX',
    min='MIN',
)


class AtomTypeEState(EStateBase):
    r'''
    atom type e-state descriptor

    Parameters:
        type(str):
            * 'count'
            * 'sum'
            * 'max'
            * 'min'

        estate(str): e-state atom type

    Returns:
        int or float: e-state value
    '''

    @classmethod
    def preset(cls):
        return (
            cls(a, t)
            for a in ['count', 'sum', 'max', 'min']
            for t in es_types
        )

    def __str__(self):
        aggr = aggr_name_dict[self.type]

        return aggr + self.estate

    descriptor_keys = 'type', 'estate'

    def __init__(self, type='count', estate='sLi'):
        assert type in ['count', 'sum', 'max', 'min']
        assert estate in es_type_set

        self.type = type
        self.estate = estate

    @property
    def dependencies(self):
        return dict(E=EStateCache())

    def calculate(self, mol, E):
        if self.type == 'count':
            return reduce(lambda a, b: a + b, E[0]).count(self.estate)

        indices = map(
            lambda e: e[1],
            filter(lambda e: self.estate in e[0], zip(*E))
        )

        try:
            return float(getattr(builtins, self.type)(indices))
        except ValueError:  # min, max to empty list
            return nan
