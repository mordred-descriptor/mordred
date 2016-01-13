from .._base import Descriptor
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


class AtomTypeEState(EStateBase):
    @classmethod
    def preset(cls):
        return (
            cls(a, t)
            for a in ['count', 'sum', 'max', 'min']
            for t in es_types
        )

    @property
    def dependencies(self):
        return dict(
            E=EStateCache.make_key()
        )

    @property
    def descriptor_name(self):
        if self.aggrigate == 'count':
            aggr = 'n'
        else:
            aggr = self.aggrigate

        return aggr + self.atom_type[:1].upper() + self.atom_type[1:]

    def __init__(self, aggrigate='count', atom_type='sLi'):
        assert aggrigate in ['count', 'sum', 'max', 'min']
        assert atom_type in es_type_set

        self.aggrigate = aggrigate
        self.atom_type = atom_type

    @property
    def descriptor_key(self):
        return self.make_key(
            self.aggrigate,
            self.atom_type,
        )

    def calculate(self, mol, E):
        if self.aggrigate == 'count':
            return reduce(lambda a, b: a + b, E[0]).count(self.atom_type)

        indices = map(
            lambda e: e[1],
            filter(lambda e: self.atom_type in e[0], zip(*E))
        )

        try:
            return float(getattr(builtins, self.aggrigate)(indices))
        except ValueError:  # min, max to empty list
            return nan

_descriptors = [AtomTypeEState]
__all__ = [d.__name__ for d in _descriptors]
