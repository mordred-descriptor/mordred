from enum import IntEnum
from functools import reduce

from rdkit.Chem import EState

from ._base import Descriptor
from ._util import parse_enum

try:
    import builtins
except ImportError:
    import __builtin__ as builtins


__all__ = ("AtomTypeEState",)


es_types = (
    "sLi",
    "ssBe",
    "ssssBe",
    "ssBH",
    "sssB",
    "ssssB",
    "sCH3",
    "dCH2",
    "ssCH2",
    "tCH",
    "dsCH",
    "aaCH",
    "sssCH",
    "ddC",
    "tsC",
    "dssC",
    "aasC",
    "aaaC",
    "ssssC",
    "sNH3",
    "sNH2",
    "ssNH2",
    "dNH",
    "ssNH",
    "aaNH",
    "tN",
    "sssNH",
    "dsN",
    "aaN",
    "sssN",
    "ddsN",
    "aasN",
    "ssssN",
    "sOH",
    "dO",
    "ssO",
    "aaO",
    "sF",
    "sSiH3",
    "ssSiH2",
    "sssSiH",
    "ssssSi",
    "sPH2",
    "ssPH",
    "sssP",
    "dsssP",
    "sssssP",
    "sSH",
    "dS",
    "ssS",
    "aaS",
    "dssS",
    "ddssS",
    "sCl",
    "sGeH3",
    "ssGeH2",
    "sssGeH",
    "ssssGe",
    "sAsH2",
    "ssAsH",
    "sssAs",
    "sssdAs",
    "sssssAs",
    "sSeH",
    "dSe",
    "ssSe",
    "aaSe",
    "dssSe",
    "ddssSe",
    "sBr",
    "sSnH3",
    "ssSnH2",
    "sssSnH",
    "ssssSn",
    "sI",
    "sPbH3",
    "ssPbH2",
    "sssPbH",
    "ssssPb",
)

es_type_set = set(es_types)


class EStateBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False


class EStateCache(EStateBase):
    __slots__ = ()

    def parameters(self):
        return ()

    def calculate(self):
        return EState.TypeAtoms(self.mol), EState.EStateIndices(self.mol)


class AggrType(IntEnum):
    __slots__ = ()

    count = 1
    sum = 2
    max = 3
    min = 4

    @property
    def as_argument(self):
        return self.name

    def description(self):
        d = {self.count: "number", self.sum: "sum", self.max: "max", self.min: "min"}
        return d[self]


aggr_names = (
    (AggrType.count, "N"),
    (AggrType.sum, "S"),
    (AggrType.max, "MAX"),
    (AggrType.min, "MIN"),
)


aggr_name_dict = dict(aggr_names)


class AtomTypeEState(EStateBase):
    r"""atom type e-state descriptor.

    :type type: str
    :param type: one of aggr_types

    :type estate: str
    :param estate: one of es_types

    :returns: NaN when type in ["min", "max"] and :math:`N_{\rm X} = 0`

    References
        * :doi:`10.1021/ci00028a014`

    """

    since = "1.0.0"
    __slots__ = ("_type", "_estate")

    def description(self):
        return "{} of {}".format(self._type.description(), self._estate)

    aggr_types = tuple(a.name for a in AggrType)

    es_types = es_types

    @classmethod
    def preset(cls, version):
        return (cls(a, t) for a in AggrType for t in es_types)

    def __str__(self):
        aggr = aggr_name_dict[self._type]

        return aggr + self._estate

    def parameters(self):
        return self._type, self._estate

    def __init__(self, type="count", estate="sLi"):
        assert estate in es_type_set

        self._type = parse_enum(AggrType, type)
        self._estate = estate

    def dependencies(self):
        return {"E": EStateCache()}

    def calculate(self, E):
        if self._type == AggrType.count:
            return reduce(lambda a, b: a + b, E[0]).count(self._estate)

        indices = map(lambda e: e[1], filter(lambda e: self._estate in e[0], zip(*E)))

        with self.rethrow_na(ValueError):
            return getattr(builtins, self._type.name)(indices)

    @property
    def rtype(self):
        r"""Return type.

        * "count": :py:class:`int`
        * "other": :py:class:`float`
        """
        return int if self._type == AggrType.count else float

    _extra_docs = "aggr_types", "es_types"
