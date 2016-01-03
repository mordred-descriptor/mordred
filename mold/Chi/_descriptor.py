from .._base import Descriptor
from .. import _atomic_property
from rdkit import Chem
from networkx import Graph
from collections import namedtuple
from enum import Enum


class ChiType(Enum):
    path = 1
    cluster = 2
    path_cluster = 3
    chain = 4


def _parse_chi_type(a):
    if isinstance(a, str):
        return getattr(ChiType, a, None)
    else:
        return ChiType(a)


class _dfs(object):
    def __init__(self, G):
        self.G = G
        self.visited = set()
        self.vis_edges = set()
        self.is_chain = False
        self.degrees = set()

    @classmethod
    def _edge_key(self, u, v):
        return min(u, v), max(u, v)

    def _dfs(self, u):
        self.visited.add(u)
        self.degrees.add(self.G.degree(u))

        for v in self.G.neighbors_iter(u):
            ek = self._edge_key(u, v)
            if v not in self.visited:
                self.vis_edges.add(ek)
                self._dfs(v)
            elif ek not in self.vis_edges:
                self.vis_edges.add(ek)
                self.is_chain = True

    def __call__(self):
        self._dfs(next(self.G.nodes_iter()))

        if self.is_chain:
            return ChiType.chain
        elif not self.degrees - set([1, 2]):
            return ChiType.path
        elif 2 in self.degrees:
            return ChiType.path_cluster
        else:
            return ChiType.cluster


class ChiBase(Descriptor):
    explicit_hydrogens = False


ChiBonds = namedtuple('ChiBonds', 'chain path path_cluster cluster')


class ChiCache(ChiBase):
    @property
    def descriptor_key(self):
        return self.make_key(self.length)

    @property
    def dependencies(self):
        return {}

    def __init__(self, length):
        self.length = length

    def calculate(self, mol):
        chain = list()
        path = list()
        path_cluster = list()
        cluster = list()
        for bonds in Chem.FindAllSubgraphsOfLengthN(mol, self.length):

            G = Graph()
            nodes = set()
            for bond in (mol.GetBondWithIdx(i) for i in bonds):
                a = bond.GetBeginAtomIdx()
                b = bond.GetEndAtomIdx()
                G.add_edge(a, b)
                nodes.add(a)
                nodes.add(b)

            typ = _dfs(G)()
            if typ == ChiType.chain:
                chain.append(nodes)
            elif typ == ChiType.path:
                path.append(nodes)
            elif typ == ChiType.path_cluster:
                path_cluster.append(nodes)
            else:
                cluster.append(nodes)

        return ChiBonds(chain, path, path_cluster, cluster)

_chi_type_dict = {
    ChiType.path: 'P',
    ChiType.chain: 'CH',
    ChiType.path_cluster: 'PC',
    ChiType.cluster: 'C'
}

_attr_dict = dict(σ='S', σv='V')


_sigmas = ['σ', 'σv']


class Chi(ChiBase):
    descriptor_defaults =\
        [(ChiType.chain, l, a) for a in _sigmas for l in range(3, 8)] +\
        [(ChiType.cluster, l, a) for a in _sigmas for l in range(3, 7)] +\
        [(ChiType.path_cluster, l, a) for a in _sigmas for l in range(4, 7)] +\
        [(ChiType.path, l, a, m) for a in _sigmas for m in [False, True] for l in range(8)]

    @property
    def descriptor_name(self):
        attr = _attr_dict.get(self.attr_name, self.attr_name)
        ct = _chi_type_dict[self.chi_type]
        p = 'A' if self.averaged else ''

        return '{}{}{}-{}'.format(p, attr, ct, self.length)

    @property
    def descriptor_key(self):
        return self.make_key(self.chi_type, self.length, self.attribute, self.averaged)

    @property
    def dependencies(self):
        chi = ChiCache.make_key(self.length) if self.length > 0 else None
        return dict(chi=chi)

    def __init__(self, chi_type=ChiType.path, length=0, attribute='σ', averaged=False):
        self.length = length
        self.attr_name, self.attribute = _atomic_property.getter(attribute)
        self.chi_type = _parse_chi_type(chi_type)
        self.averaged = averaged

    def calculate(self, mol, chi):
        if self.length <= 0:
            chi = ChiBonds([], [{a.GetIdx()} for a in mol.GetAtoms()], [], [])

        x = 0
        node_sets = getattr(chi, self.chi_type.name)
        for nodes in node_sets:
            c = 1
            for node in nodes:
                c *= self.attribute(mol.GetAtomWithIdx(node))

            x += c ** (-0.5)

        if self.averaged:
            x /= len(node_sets) or 1

        return x


_descriptors = [Chi]
__all__ = [d.__name__ for d in _descriptors]
