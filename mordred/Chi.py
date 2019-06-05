from enum import IntEnum
from itertools import chain
from collections import namedtuple, defaultdict

from rdkit import Chem

from ._base import Descriptor
from ._util import parse_enum
from ._atomic_property import AtomicProperty

__all__ = ("Chi",)


class ChiType(IntEnum):
    __slots__ = ()

    path = 1
    cluster = 2
    path_cluster = 3
    chain = 4

    @property
    def as_argument(self):
        return self.name

    @property
    def short(self):
        _short_dict = {
            self.path: "p",
            self.chain: "ch",
            self.path_cluster: "pc",
            self.cluster: "c",
        }

        return _short_dict[self]

    @property
    def long(self):
        _long_dict = {
            self.path: "Chi path",
            self.chain: "Chi chain",
            self.path_cluster: "Chi path-cluster",
            self.cluster: "Chi cluster",
        }

        return _long_dict[self]


class DFS(object):
    __slots__ = (
        "mol",
        "visited",
        "vis_edges",
        "is_chain",
        "degrees",
        "bonds",
        "neighbors",
    )

    def __init__(self, mol):
        self.mol = mol
        self.visited = set()
        self.vis_edges = set()
        self.degrees = set()
        self.neighbors = defaultdict(set)

        self.bonds = [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in self.mol.GetBonds()
        ]

    def reset(self, use_bonds):
        ns = self.neighbors
        bs = self.bonds

        self.is_chain = False
        self.visited.clear()
        self.vis_edges.clear()
        self.degrees.clear()
        ns.clear()

        for i in range(len(use_bonds)):
            a, b = bs[use_bonds[i]]
            ns[a].add(b)
            ns[b].add(a)

        self.neighbors = ns

    @property
    def nodes(self):
        return list(self.neighbors.keys())

    def _dfs(self, u):
        neighbors = self.neighbors[u]
        self.visited.add(u)
        self.degrees.add(len(neighbors))

        for v in neighbors:
            ek = (v, u) if u > v else (u, v)

            if v not in self.visited:
                self.vis_edges.add(ek)
                self._dfs(v)

            elif ek not in self.vis_edges:
                self.vis_edges.add(ek)
                self.is_chain = True

    def __call__(self):
        self._dfs(next(iter(self.neighbors.keys())))

        if self.is_chain:
            t = ChiType.chain
        elif not self.degrees - {1, 2}:
            t = ChiType.path
        elif 2 in self.degrees:
            t = ChiType.path_cluster
        else:
            t = ChiType.cluster

        return t


class ChiBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False


ChiBonds = namedtuple("ChiBonds", "chain path path_cluster cluster")


class ChiCache(ChiBase):
    __slots__ = ("_order",)

    def parameters(self):
        return (self._order,)

    def __init__(self, order):
        self._order = order

    def calculate(self):
        chain = []
        path = []
        path_cluster = []
        cluster = []

        dfs = DFS(self.mol)
        for bonds in Chem.FindAllSubgraphsOfLengthN(self.mol, self._order):
            dfs.reset(bonds)
            typ = dfs()
            nodes = dfs.nodes

            if typ == ChiType.chain:
                chain.append(nodes)
            elif typ == ChiType.path:
                path.append(nodes)
            elif typ == ChiType.path_cluster:
                path_cluster.append(nodes)
            else:
                cluster.append(nodes)

        return ChiBonds(chain, path, path_cluster, cluster)


class Chi(ChiBase):
    r"""chi descriptor.

    :type type: str
    :param type: one of chi_types

    :type prop: str or function
    :param prop: :ref:`atomic_properties`

    :type averaged: bool
    :param averaged: averaged by number of subgraphs

    :returns: NaN when

        * any atomic properties <= 0
        * averaged and :math:`N_{\chi} = 0`
    """

    since = "1.0.0"
    __slots__ = ("_type", "_order", "_prop", "_averaged")

    chi_types = tuple(t.name for t in ChiType)

    _deltas = ["d", "dv"]

    def description(self):
        return "{}-ordered {}{} weighted by {}".format(
            self._order,
            "averaged " if self._averaged else "",
            self._type.long,
            self._prop.get_long(),
        )

    @classmethod
    def preset(cls, version):
        return chain(
            (cls(ChiType.chain, l, a) for a in cls._deltas for l in range(3, 8)),
            (cls(ChiType.cluster, l, a) for a in cls._deltas for l in range(3, 7)),
            (cls(ChiType.path_cluster, l, a) for a in cls._deltas for l in range(4, 7)),
            (
                cls(ChiType.path, l, a, m)
                for a in cls._deltas
                for m in [False, True]
                for l in range(8)
            ),
        )

    def __str__(self):
        prop = self._prop.as_argument
        ct = self._type.short
        averaged = "A" if self._averaged else ""

        return "{}X{}-{}{}".format(averaged, ct, self._order, prop)

    def parameters(self):
        return self._type, self._order, self._prop, self._averaged

    def __init__(self, type="path", order=0, prop="d", averaged=False):
        self._type = parse_enum(ChiType, type)
        self._order = order
        self._prop = AtomicProperty(self.explicit_hydrogens, prop)
        self._averaged = averaged

    def dependencies(self):
        d = {"P": self._prop}
        if self._order > 0:
            d["chi"] = ChiCache(self._order)

        return d

    def calculate(self, P, chi=None):
        if self._order <= 0:
            chi = ChiBonds([], [{a.GetIdx()} for a in self.mol.GetAtoms()], [], [])

        x = 0.0
        node_sets = getattr(chi, self._type.name)
        for nodes in node_sets:
            c = 1
            for node in nodes:
                c *= P[node]

            if c <= 0:
                self.fail(ValueError("some properties less then or equal to 0"))

            x += c ** -0.5

        if self._averaged:
            with self.rethrow_zerodiv():
                x /= len(node_sets)

        return x

    rtype = float

    _extra_docs = ("chi_types",)
