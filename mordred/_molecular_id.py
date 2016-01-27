import math

from networkx import Graph

from rdkit import Chem

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


table = Chem.GetPeriodicTable()


class MolecularIdBase(Descriptor):
    explicit_hydrogens = False
    require_connected = True


class AtomicIds(MolecularIdBase):
    __slots__ = ()

    def calculate(self, mol):
        aid = AtomicId(mol)
        return [
            1 + aid.get_atomic_id(i) / 2.0
            for i in range(mol.GetNumAtoms())
        ]


class MolecularId(MolecularIdBase):
    r"""molecular id descriptor.

    :type type: :py:class:`str` or :py:class`int`
    :param type: target of atomic id source

        * 'any': normal molecular id(sum of all atomic id)
        * 'X': sum of halogen atomic id
        * str: atomic symbol
        * int: atomic number

    :type averaged: bool
    :param averaged: averaged by number of atoms

    :rtype: float
    """

    @classmethod
    def preset(cls):
        return (
            cls(s, a)
            for s in ['any', 'hetero', 'C', 'N', 'O', 'X']
            for a in [False, True]
        )

    def __str__(self):
        n = 'AMID' if self._averaged else 'MID'
        if self._type != 'any':
            n = '{}_{}'.format(n, self._type)

        return n

    __slots__ = ('_orig_type', '_averaged',)

    def __init__(self, type='any', averaged=False):
        self._orig_type = self._type = type
        self._averaged = averaged

        if isinstance(type, str) and type not in ['any', 'hetero', 'X']:
            type = table.GetAtomicNumber(type)

        if type == 'any':
            self._check = lambda _: True
        elif type == 'hetero':
            self._type = 'h'
            self._check = lambda a: a not in set([1, 6])
        elif self._type == 'X':
            self._check = lambda a: a in set([9, 17, 35, 53, 85, 117])
        else:
            self._check = lambda a: a == type

    def dependencies(self):
        return dict(aids=AtomicIds())

    def calculate(self, mol, aids):
        v = sum(
            aid
            for aid, atom in zip(aids, mol.GetAtoms())
            if self._check(atom.GetAtomicNum())
        )

        if self._averaged:
            v /= mol.GetNumAtoms()

        return v
