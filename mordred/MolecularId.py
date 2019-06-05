import math

from six import integer_types
from networkx import Graph

from ._base import Descriptor
from ._atomic_property import GetAtomicNumber, GetElementSymbol, halogen

__all__ = ("MolecularId",)


class AtomicId(object):
    __slots__ = ("G", "lim", "start", "id", "visited", "weights")

    def __init__(self, mol, eps):
        G = Graph()

        G.add_nodes_from(a.GetIdx() for a in mol.GetAtoms())

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
            w = d["weight"] * self.weights[-1]
            self.weights.append(w)

            self.id += 1.0 / math.sqrt(w)
            if w < self.lim:
                self._search(v)

            self.visited.remove(v)
            self.weights.pop()

    def __call__(self):
        return [self.get_atomic_id(i) for i in range(self.G.number_of_nodes())]


class MolecularIdBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False
    require_connected = True

    def parameters(self):
        return (self._eps,)


class AtomicIds(MolecularIdBase):
    __slots__ = ("_eps",)

    def __init__(self, eps=1e-10):
        self._eps = eps

    def calculate(self):
        aid = AtomicId(self.mol, self._eps)
        return [1 + aid.get_atomic_id(i) / 2.0 for i in range(self.mol.GetNumAtoms())]


class MolecularId(MolecularIdBase):
    r"""molecular id descriptor.

    :type type: :py:class:`str` or :py:class:`int`
    :param type: target of atomic id source

        * "any": normal molecular id(sum of all atomic id)
        * "X": sum of halogen atomic id
        * str: atomic symbol
        * int: atomic number

    :type averaged: bool
    :param averaged: averaged by number of atoms

    :type _eps: float
    :param _eps: internally used
    """

    since = "1.0.0"
    __slots__ = ("_orig_type", "_averaged", "_eps", "_type", "_check")

    def description(self):
        if self._type == "any":
            t = ""
        elif self._type == "X":
            t = " on halogen atoms"
        else:
            e = self._type
            if isinstance(e, integer_types):
                e = GetElementSymbol(e)

            t = " on {} atoms".format(e)

        return "{}molecular ID{}".format("averaged " if self._averaged else "", t)

    @classmethod
    def preset(cls, version):
        return (
            cls(s, a)
            for s in ["any", "hetero", "C", "N", "O", "X"]
            for a in [False, True]
        )

    def __str__(self):
        n = "AMID" if self._averaged else "MID"
        if self._type != "any":
            n = "{}_{}".format(n, self._type)

        return n

    def parameters(self):
        return self._orig_type, self._averaged, self._eps

    def __init__(self, type="any", averaged=False, _eps=1e-10):
        self._orig_type = self._type = type
        self._averaged = averaged
        self._eps = _eps

        if isinstance(type, str) and type not in ["any", "hetero", "X"]:
            type = GetAtomicNumber(type)

        if type == "any":
            self._check = lambda _: True
        elif type == "hetero":
            self._type = "h"
            self._check = lambda a: a not in {1, 6}
        elif self._type == "X":
            self._check = lambda a: a in halogen
        else:
            self._check = lambda a: a == type

    def dependencies(self):
        return {"aids": AtomicIds(self._eps)}

    def calculate(self, aids):
        v = float(
            sum(
                aid
                for aid, atom in zip(aids, self.mol.GetAtoms())
                if self._check(atom.GetAtomicNum())
            )
        )

        if self._averaged:
            v /= self.mol.GetNumAtoms()

        return v

    rtype = float
