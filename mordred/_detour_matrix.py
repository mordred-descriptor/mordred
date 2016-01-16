import networkx
import numpy as np
from ._base import Descriptor
from . import _matrix_attributes as ma


class LongestSimplePath(object):
    def __init__(self, G, weight=None):
        self.G = G
        self.N = G.number_of_nodes()
        self.neighbors = {n: [(v, d.get(weight, 1.0))
                              for (v, d) in G[n].items()]
                          for n in G.nodes()}

    def _start(self, s):
        self.start = s
        self.result = {n: 0 for n in self.G.nodes_iter()}
        self.visited = set()
        self.distance = 0.0
        self._search(s)
        return self.result

    def _search(self, u):
        self.visited.add(u)
        for v, w in self.neighbors[u]:
            if v in self.visited:
                continue

            self.visited.add(v)
            self.distance += w

            d = self.distance
            if d > self.result[v]:
                self.result[v] = d

            if v != self.start:
                self._search(v)

            self.visited.remove(v)
            self.distance -= w

    def __call__(self):
        return {(min(s, g), max(s, g)): w
                for s in self.G.nodes_iter()
                for g, w in self._start(s).items()}


class calc_detour(object):
    def __init__(self, G, weight='weight'):
        self.N = G.number_of_nodes()
        self.Q = []
        for bcc in networkx.biconnected_component_subgraphs(G, False):
            lsp = LongestSimplePath(bcc, weight)()
            nodes = set()
            for a, b in lsp:
                nodes.add(a)
                nodes.add(b)
            self.Q.append((nodes, lsp))

        self.nodes, self.C = self.Q.pop()

    def merge(self):
        for i in range(1, len(self.Q) + 1):
            ns, lsp = self.Q[-i]
            common = ns & self.nodes
            if len(common) == 0:
                continue
            elif len(common) > 1:
                raise ValueError('bug: multiple common nodes.')

            common = common.pop()
            self.Q.pop(-i)
            for n in ns:
                self.nodes.add(n)
            break

        def calc_weight(i, j):
            ij = (i, j)
            if ij in self.C:
                return self.C[ij]
            elif ij in lsp:
                return lsp[ij]
            elif i == j == common:
                return max(lsp[ij], self.C[ij])

            ic = (min(i, common), max(i, common))
            jc = (min(j, common), max(j, common))

            if ic in self.C and jc in lsp:
                return self.C[ic] + lsp[jc]
            elif jc in self.C and ic in lsp:
                return self.C[jc] + lsp[ic]
            else:
                raise ValueError('bug: unknown weight')

        self.C = {(i, j): calc_weight(i, j)
                  for i in self.nodes
                  for j in self.nodes
                  if i <= j}

    def __call__(self):
        while self.Q:
            self.merge()

        result = np.empty((self.N, self.N))
        for i, j in ((i, j) for i in range(self.N) for j in range(i, self.N)):
            result[i, j] = self.C[(i, j)]
            result[j, i] = self.C[(i, j)]

        return result


class DetourMatrixBase(Descriptor):
    explicit_hydrogens = False


class detour_matrix(DetourMatrixBase):
    def calculate(self, mol):
        G = networkx.Graph()
        G.add_edges_from(
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            for b in mol.GetBonds()
        )

        return calc_detour(G)()


class DetourMatrix(DetourMatrixBase):
    r'''
    detour matrix descriptor

    Parameters:
        method(str): matrix aggregate method

    Returns:
        float: result
    '''

    @classmethod
    def preset(cls):
        return map(cls, ma.methods)

    def __str__(self):
        return '{}_Dt'.format(self.method.__name__)

    descriptor_keys = 'method',

    def __init__(self, method='SpMax'):
        self.method = ma.get_method(method)

    @property
    def dependencies(self):
        return dict(
            result=self.method(
                detour_matrix(),
                self.explicit_hydrogens,
                self.gasteiger_charges,
                self.kekulize,
            )
        )

    def calculate(self, mol, result):
        return result
