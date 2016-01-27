import math

from networkx import Graph

from ._base import Descriptor


class AtomicId(object):
    def __init__(self, mol, eps=1e-10):
        G = Graph()

        for bond in mol.GetBonds():
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()

            w = a.GetDegree() * b.GetDegree()

            G.add_edge(a.GetIdx(), b.GetIdx(), weight=w)

        self.G = G
        self.lim = int(1.0 / (eps ** 2))

    def get_atomic_id(self, s):
        self.start = s
        self.id = 0.0
        self.visited = set()
        self.weights = [1]
        self._search(s)
        return self.id

    def _search(self, u):
        self.visited.add(u)

        for v, d in self.G[u].items():
            if v in self.visited:
                continue

            self.visited.add(v)
            w = d['weight'] * self.weights[-1]
            self.weights.append(w)

            self.id += 1.0 / math.sqrt(w)
            if w < self.lim:
                self._search(v)

            self.visited.remove(v)
            self.weights.pop()

    def __call__(self):
        return [
            self.get_atomic_id(i)
            for i in range(self.G.number_of_nodes())
        ]


class AtomicIds(Descriptor):
    __slots__ = ()

    explicit_hydrogens = False

    def calculate(self, mol):
        aid = AtomicId(mol)
        return [
            1 + aid.get_atomic_id(i) / 2.0
            for i in range(mol.GetNumAtoms())
        ]


class MolecularId(Descriptor):
    r"""molecular id descriptor.

    :type averaged: bool
    :param averaged: averaged by number of atoms

    :rtype: float
    """

    @classmethod
    def preset(cls):
        return (cls(b) for b in [False, True])

    explicit_hydrogens = False

    def __str__(self):
        return 'AMID' if self._averaged else 'MID'

    __slots__ = ('_averaged',)

    def __init__(self, averaged=False):
        self._averaged = averaged

    def dependencies(self):
        return dict(aid=AtomicIds())

    def calculate(self, mol, aid):
        v = sum(aid)

        if self._averaged:
            v /= mol.GetNumAtoms()

        return v
