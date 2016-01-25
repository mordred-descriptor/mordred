from enum import IntEnum
from functools import reduce

from numpy import nan

from rdkit.Chem import EState

from ._base import Descriptor, parse_enum

try:
    import builtins
except ImportError:
    import __builtin__ as builtins


es_types = (
    'sLi', 'ssBe', 'ssssBe', 'ssBH', 'sssB', 'ssssB', 'sCH3', 'dCH2', 'ssCH2',
    'tCH', 'dsCH', 'aaCH', 'sssCH', 'ddC', 'tsC', 'dssC', 'aasC', 'aaaC',
    'ssssC', 'sNH3', 'sNH2', 'ssNH2', 'dNH', 'ssNH', 'aaNH', 'tN', 'sssNH',
    'dsN', 'aaN', 'sssN', 'ddsN', 'aasN', 'ssssN', 'sOH', 'dO', 'ssO', 'aaO',
    'sF', 'sSiH3', 'ssSiH2', 'sssSiH', 'ssssSi', 'sPH2', 'ssPH', 'sssP',
    'dsssP', 'sssssP', 'sSH', 'dS', 'ssS', 'aaS', 'dssS', 'ddssS', 'sCl',
    'sGeH3', 'ssGeH2', 'sssGeH', 'ssssGe', 'sAsH2', 'ssAsH', 'sssAs', 'sssdAs',
    'sssssAs', 'sSeH', 'dSe', 'ssSe', 'aaSe', 'dssSe', 'ddssSe', 'sBr', 'sSnH3',
    'ssSnH2', 'sssSnH', 'ssssSn', 'sI', 'sPbH3', 'ssPbH2', 'sssPbH', 'ssssPb'
)

es_type_set = set(es_types)


class EStateBase(Descriptor):
    explicit_hydrogens = False


class EStateCache(EStateBase):
    __slots__ = ()

    def calculate(self, mol):
        return EState.TypeAtoms(mol), EState.EStateIndices(mol)


class AggrType(IntEnum):
    count = 1
    sum = 2
    max = 3
    min = 4

aggr_names = (
    (AggrType.count, 'N'),
    (AggrType.sum, 'S'),
    (AggrType.max, 'MAX'),
    (AggrType.min, 'MIN'),
)


aggr_name_dict = dict(aggr_names)


class AtomTypeEState(EStateBase):
    r"""atom type e-state descriptor.

    :type type: str
    :param type: one of aggr_types

    :type estate: str
    :param estate: one of es_types

    :rtype: int('count') or float(other)
    :returns: NaN when type in ['min', 'max'] and :math:`N_{\rm X} = 0`

    References
        * :cite:`10.1021/ci00028a014`
    """

    aggr_types = tuple(a.name for a in AggrType)

    es_types = es_types

    @classmethod
    def preset(cls):
        return (
            cls(a, t)
            for a in AggrType
            for t in es_types
        )

    def __str__(self):
        aggr = aggr_name_dict[self._type]

        return aggr + self._estate

    __slots__ = ('_type', '_estate',)

    def __init__(self, type='count', estate='sLi'):
        assert estate in es_type_set

        self._type = parse_enum(AggrType, type)
        self._estate = estate

    def dependencies(self):
        return dict(E=EStateCache())

    def calculate(self, mol, E):
        if self._type == AggrType.count:
            return reduce(lambda a, b: a + b, E[0]).count(self._estate)

        indices = map(
            lambda e: e[1],
            filter(lambda e: self._estate in e[0], zip(*E))
        )

        try:
            return float(getattr(builtins, self._type.name)(indices))
        except ValueError:  # min, max to empty list
            return nan
