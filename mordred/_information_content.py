from itertools import groupby

import numpy as np

from ._base import Descriptor
from ._common import DistanceMatrix


class BFSTree(object):
    __slots__ = ('mol', 'tree', 'order', 'visited',)

    def __init__(self, mol, i):
        self.mol = mol
        self.tree = {i: ()}
        self.order = 0
        self.visited = set([i])

    def expand(self):
        self._expand(self.tree)
        self.order += 1

    def _expand(self, tree):
        for src, dst in list(tree.items()):
            self.visited.add(src)

            if dst is ():
                tree[src] = {
                    n.GetIdx(): ()
                    for n in self.mol.GetAtomWithIdx(src).GetNeighbors()
                    if n.GetIdx() not in self.visited
                }

            else:
                self._expand(dst)

    @property
    def code(self):
        return tuple(sorted(self._code(self.tree, None, ())))

    def _code(self, tree, before, trail):
        if len(tree) == 0:
            yield trail

        else:
            for src, dst in tree.items():

                code = []
                if before is not None:
                    code.append(self.mol.GetBondBetweenAtoms(before, src).GetBondType())

                a = self.mol.GetAtomWithIdx(src)
                code.append(a.GetAtomicNum())
                code.append(a.GetDegree())

                nxt = tuple(list(trail) + code)
                for t in self._code(dst, src, nxt):
                    yield t


def neighborhood_code(mol, i, order):
    if order == 0:
        return mol.GetAtomWithIdx(i).GetAtomicNum()

    tree = BFSTree(mol, i)
    for _ in range(order):
        tree.expand()
    return tree.code


class InformationContentBase(Descriptor):
    kekulize = True
    require_connected = False


class Ag(InformationContentBase):
    __slots__ = ('_order',)

    def __init__(self, order):
        self._order = order

    def dependencies(self):
        return dict(
            D=DistanceMatrix(self.explicit_hydrogens)
        )

    def calculate(self, mol, D):
        atoms = [neighborhood_code(mol, i, self._order) for i in range(mol.GetNumAtoms())]
        ad = {a: i for i, a in enumerate(atoms)}
        Ags = [(k, sum(1 for _ in g)) for k, g in groupby(sorted(atoms))]
        return\
            np.array([ad[k] for k, _ in Ags]),\
            np.array([ag for _, ag in Ags])


class InformationContent(InformationContentBase):
    r"""information content descriptor.

    :type type: str
    :param type: one of ic_types

    :type order: int
    :param order: order(number of edge) of subgraph

    :rtype: float
    """

    ic_types = ('', 'T', 'S', 'C', 'B', 'M', 'ZM')

    @classmethod
    def preset(cls):
        return (
            cls(t, o)
            for t in cls.ic_types
            for o in range(6)
        )

    def __str__(self):
        return '{}IC{}'.format(self._type, self._order)

    __slots__ = ('_type', '_order',)

    def __init__(self, type='', order=0):
        assert type in self.ic_types
        self._type = type
        self._order = order

    def dependencies(self):
        if self._type in ['T', 'S', 'C', 'B']:
            return dict(
                ICm=self.__class__('', self._order)
            )
        else:
            return dict(
                iAgs=Ag(self._order)
            )

    def calculate(self, mol, iAgs=None, ICm=None):
        N = mol.GetNumAtoms()
        if self._type in ['', 'M', 'ZM']:
            ids = iAgs[0]
            Ags = iAgs[1]

            if self._type == '':
                w = 1
            elif self._type == 'M':
                w = np.vectorize(lambda i: mol.GetAtomWithIdx(int(i)).GetMass())(ids)
            else:
                w = Ags * np.vectorize(lambda i: mol.GetAtomWithIdx(int(i)).GetAtomicNum())(ids)

            e = np.vectorize(lambda Ag: entropy_term(float(Ag) / N))(Ags)
            return -np.sum(w * e)

        if self._type == 'T':
            return N * ICm
        elif self._type == 'S':
            return ICm / np.log2(N)
        elif self._type == 'C':
            return np.log2(N) - ICm
        elif self._type == 'B':
            bts = sum(b.GetBondTypeAsDouble() for b in mol.GetBonds())
            if bts <= 1:
                return np.nan

            return ICm / np.log2(bts)


def entropy_term(a):
    return a * np.log2(a)
