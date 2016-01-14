from .._base import Descriptor
from ..Ring._descriptor import Rings
import networkx as nx


class FrameworkCache(Descriptor):
    @property
    def dependencies(self):
        return dict(
            Rs=Rings.make_key()
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
        for Ri, Rj in ((i, j) for i in range(NR) for j in range(i+1, NR)):
            Ra, Rb = R[Ri], R[Rj]
            linkers.update(i for t, i in nx.shortest_path(G, Ra, Rb) if t == 'A')

        return linkers, Rs


class Framework(Descriptor):
    r'''
    molecular framework ratio descriptor

    .. math::

        f_{\rm MF} = \frac{N_{\rm MF}}{N}

    where
    :math:`N_{\rm MF}` is number of atoms in molecular framework,
    :math:`N` is number of all atoms.

    Returns:
        float: fMF value
    '''

    descriptor_name = 'fMF'

    @property
    def dependencies(self):
        return dict(
            F=FrameworkCache.make_key()
        )

    def calculate(self, mol, F):
        linkers, rings = F
        Nmf = len(linkers) + len({i for ring in rings for i in ring})
        N = mol.GetNumAtoms()

        return float(Nmf) / float(N)

_descriptors = [Framework]
__all__ = [d.__name__ for d in _descriptors]
