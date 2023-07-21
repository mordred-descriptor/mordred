import time

import numpy as np
import networkx

from . import _matrix_attributes as ma
from ._base import Descriptor
from .error import Timeout

__all__ = ("DetourMatrix", "DetourIndex")


class LongestSimplePath(object):
    __slots__ = (
        "G",
        "N",
        "neighbors",
        "start",
        "result",
        "visited",
        "distance",
        "timeout_at",
    )

    def __init__(self, G, weight=None, timeout_at=None):
        self.G = G
        self.N = G.number_of_nodes()
        self.timeout_at = timeout_at
        self.neighbors = {
            n: [(v, d.get(weight, 1.0)) for v, d in G[n].items()] for n in G.nodes()
        }

    def _start(self, s):
        self.start = s
        self.result = {n: 0 for n in self.G.nodes()}
        self.visited = set()
        self.distance = 0.0
        self._search(s)
        return self.result

    def _search(self, u):
        if self.timeout_at < time.time():
            raise Timeout()

        self.visited.add(u)
        for v, w in self.neighbors[u]:
            if v in self.visited:
                continue

            self.visited.add(v)
            self.distance += w

            d = self.distance
            if d > self.result[v]:
                self.result[v] = d

            self._search(v)

            self.visited.remove(v)
            self.distance -= w

    def __call__(self):
        return {
            (min(s, g), max(s, g)): w
            for s in self.G.nodes()
            for g, w in self._start(s).items()
        }


class CalcDetour(object):
    __slots__ = ("N", "G", "Q", "nodes", "C", "weight", "timeout")

    def __init__(self, G, weight="weight", timeout=None):
        self.G = G
        self.N = G.number_of_nodes()
        self.Q = []
        self.weight = weight
        self.timeout = timeout

    def merge(self):
        for i in range(1, len(self.Q) + 1):
            ns, lsp = self.Q[-i]
            common = ns & self.nodes
            if len(common) == 0:
                continue
            elif len(common) > 1:
                raise ValueError("bug: multiple common nodes.")

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
                raise ValueError("bug: unknown weight")

        self.C = {
            (i, j): calc_weight(i, j) for i in self.nodes for j in self.nodes if i <= j
        }

    def __call__(self):
        timeout_at = None if self.timeout is None else time.time() + self.timeout

        for bcc in (
            self.G.subgraph(c) for c in networkx.biconnected_components(self.G)
        ):
            lsp = LongestSimplePath(bcc, self.weight, timeout_at)()
            nodes = set()
            for a, b in lsp:
                nodes.add(a)
                nodes.add(b)
            self.Q.append((nodes, lsp))

        if self.N == 1:
            return np.array([[0]])

        self.nodes, self.C = self.Q.pop()

        while self.Q:
            self.merge()

        result = np.empty((self.N, self.N))
        for i, j in ((i, j) for i in range(self.N) for j in range(i, self.N)):
            result[i, j] = self.C[(i, j)]
            result[j, i] = self.C[(i, j)]

        return result


class DetourMatrixBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False
    require_connected = True


class DetourMatrixCache(DetourMatrixBase):
    __slots__ = ()

    hermitian = True

    def parameters(self):
        return ()

    def calculate(self):
        G = networkx.Graph()
        G.add_nodes_from(a.GetIdx() for a in self.mol.GetAtoms())
        G.add_edges_from(
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in self.mol.GetBonds()
        )

        return CalcDetour(G, timeout=self.config.get("timeout", 60))()


class DetourMatrix(DetourMatrixBase):
    r"""detour matrix descriptor.

    :type type: str
    :param type: :ref:`matrix_aggregating_methods`
    """

    since = "1.0.0"
    __slots__ = ("_type",)

    def description(self):
        return "{} from detourn matrix".format(self._type.description())

    @classmethod
    def preset(cls, version):
        return map(cls, ma.methods)

    def __str__(self):
        return "{}_Dt".format(self._type.__name__)

    def parameters(self):
        return (self._type,)

    def __init__(self, type="SpMax"):
        self._type = ma.get_method(type)

    def dependencies(self):
        return {
            "result": self._type(
                DetourMatrixCache(), self.explicit_hydrogens, self.kekulize
            )
        }

    def calculate(self, result):
        return result

    rtype = float


class DetourIndex(DetourMatrixBase):
    r"""detour index descriptor.

    .. math::

        I_{\rm D} = \frac{1}{A} \sum^A_{i=1} \sum^A_{j=1} {\boldsymbol D}_{ij}

    where
    :math:`D` is detour matrix,
    :math:`A` is number of atoms.
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "detour index"

    def parameters(self):
        return ()

    @classmethod
    def preset(cls, version):
        yield cls()

    explicit_hydrogens = False

    def __str__(self):
        return self.__class__.__name__

    def dependencies(self):
        return {"D": DetourMatrixCache()}

    def calculate(self, D):
        return int(0.5 * D.sum())

    rtype = int
