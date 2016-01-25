import networkx as nx

from ._base import Descriptor
from ._ring_count import Rings


class FrameworkBase(Descriptor):
    pass


class FrameworkCache(FrameworkBase):
    __slots__ = ()

    def dependencies(self):
        return dict(
            Rs=Rings()
        )

    def calculate(self, mol, Rs):
        G = nx.Graph()
        Rd = {i: ('R', Ri) for Ri, R in enumerate(Rs) for i in R}
        R = list(set(Rd.values()))
        NR = len(R)

        for bond in mol.GetBonds():
            a = bond.GetBeginAtomIdx()
            b = bond.GetEndAtomIdx()

            a = Rd.get(a, ('A', a))
            b = Rd.get(b, ('A', b))

            G.add_edge(a, b)

        linkers = set()
        for Ri, Rj in ((i, j) for i in range(NR) for j in range(i + 1, NR)):
            Ra, Rb = R[Ri], R[Rj]
            try:
                linkers.update(i for t, i in nx.shortest_path(G, Ra, Rb) if t == 'A')
            except nx.NetworkXNoPath:
                pass

        return linkers, Rs


class Framework(FrameworkBase):
    r"""molecular framework ratio descriptor.

    .. math::

        f_{\rm MF} = \frac{N_{\rm MF}}{N}

    where
    :math:`N_{\rm MF}` is number of atoms in molecular framework,
    :math:`N` is number of all atoms.

    :rtype: float

    References
        * :cite:`10.1021/jm9602928`
    """

    __slots__ = ()

    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return 'fMF'

    def dependencies(self):
        return dict(
            F=FrameworkCache()
        )

    def calculate(self, mol, F):
        linkers, rings = F
        Nmf = len(linkers) + len({i for ring in rings for i in ring})
        N = mol.GetNumAtoms()

        return float(Nmf) / float(N)
