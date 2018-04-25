from __future__ import division

import networkx as nx

from ._base import Descriptor
from .RingCount import Rings

__all__ = ("Framework",)


class FrameworkCache(Descriptor):
    __slots__ = ()

    def parameters(self):
        return ()

    def dependencies(self):
        return {"Rs": Rings()}

    def calculate(self, Rs):
        G = nx.Graph()
        Rd = {i: ("R", Ri) for Ri, R in enumerate(Rs) for i in R}
        R = list(set(Rd.values()))
        NR = len(R)

        for bond in self.mol.GetBonds():
            a = bond.GetBeginAtomIdx()
            b = bond.GetEndAtomIdx()

            a = Rd.get(a, ("A", a))
            b = Rd.get(b, ("A", b))

            G.add_edge(a, b)

        linkers = set()
        for Ri, Rj in ((i, j) for i in range(NR) for j in range(i + 1, NR)):
            Ra, Rb = R[Ri], R[Rj]
            try:
                linkers.update(i for t, i in nx.shortest_path(G, Ra, Rb) if t == "A")
            except nx.NetworkXNoPath:
                pass

        return linkers, Rs


class Framework(Descriptor):
    r"""molecular framework ratio descriptor.

    .. math::

        f_{\rm MF} = \frac{N_{\rm MF}}{N}

    where
    :math:`N_{\rm MF}` is number of atoms in molecular framework,
    :math:`N` is number of all atoms.

    References
        * :doi:`10.1021/jm9602928`

    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "molecular framework ratio"

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return "fMF"

    def parameters(self):
        return ()

    def dependencies(self):
        return {"F": FrameworkCache()}

    def calculate(self, F):
        linkers, rings = F
        Nmf = len(linkers) + len({i for ring in rings for i in ring})
        N = self.mol.GetNumAtoms()

        return Nmf / N

    rtype = float
