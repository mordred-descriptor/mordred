from ._base import Descriptor
from . import _atomic_property
from rdkit import Chem
from networkx import Graph
from collections import namedtuple
from enum import Enum
from itertools import chain
import numpy as np


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
    descriptor_keys = 'order',

    def __init__(self, order):
        self.order = order

    def calculate(self, mol):
        chain = list()
        path = list()
        path_cluster = list()
        cluster = list()
        for bonds in Chem.FindAllSubgraphsOfLengthN(mol, self.order):

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

_prop_dict = dict(delta='S', delta_v='V')


_deltas = ['delta', 'delta_v']


class Chi(ChiBase):
    r'''
    chi descriptor

    Parameters:
        type(str):

            * 'path'
            * 'path-cluster'
            * 'cluster'
            * 'chain'

        prop(str, function): atomic property

        averaged(bool): averaged by number of subgraphs

    Returns:
        float: chi value
    '''

    @classmethod
    def preset(cls):
        return chain(
            (cls(ChiType.chain, l, a) for a in _deltas for l in range(3, 8)),
            (cls(ChiType.cluster, l, a) for a in _deltas for l in range(3, 7)),
            (cls(ChiType.path_cluster, l, a) for a in _deltas for l in range(4, 7)),
            (cls(ChiType.path, l, a, m) for a in _deltas for m in [False, True] for l in range(8)),
        )

    def __str__(self):
        prop = _prop_dict.get(self.prop_name, self.prop_name)
        ct = _chi_type_dict[self.type]
        p = 'A' if self.averaged else ''

        return '{}{}{}-{}'.format(p, prop, ct, self.order)

    @property
    def gasteiger_charges(self):
        return getattr(self.prop, 'gasteiger_charges', False)

    descriptor_keys = 'type', 'order', 'prop', 'averaged'

    def __init__(self, type='path', order=0, prop='delta', averaged=False):
        self.order = order
        self.prop_name, self.prop = _atomic_property.getter(prop, self.explicit_hydrogens)
        self.type = _parse_chi_type(type)
        self.averaged = averaged

    @property
    def dependencies(self):
        if self.order > 0:
            return dict(chi=ChiCache(self.order))

    def calculate(self, mol, chi=None):
        if self.order <= 0:
            chi = ChiBonds([], [{a.GetIdx()} for a in mol.GetAtoms()], [], [])

        x = 0.0
        node_sets = getattr(chi, self.type.name)
        props = [self.prop(a) for a in mol.GetAtoms()]
        for nodes in node_sets:
            c = 1
            for node in nodes:
                c *= props[node]

            x += c ** -0.5

        if self.averaged:
            x /= len(node_sets) or np.nan

        return x
