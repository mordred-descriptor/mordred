from ._base import Descriptor
from ._common import DistanceMatrix
import numpy as np
from itertools import groupby
import math


class BFSTree(object):
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


class Ag(Descriptor):
    @property
    def dependencies(self):
        return dict(D=DistanceMatrix.make_key(
            self.explicit_hydrogens,
            False,
            False,
        ))

    def __init__(self, order):
        self.order = order

    @property
    def descriptor_key(self):
        return self.make_key(self.order)

    def calculate(self, mol, D):
        atoms = [neighborhood_code(mol, i, self.order) for i in range(mol.GetNumAtoms())]
        ad = {a: i for i, a in enumerate(atoms)}
        Ags = [(k, sum(1 for _ in g)) for k, g in groupby(sorted(atoms))]
        return\
            np.array([ad[k] for k, _ in Ags]),\
            np.array([ag for _, ag in Ags])


class InformationContent(Descriptor):
    r'''
    information content descriptor

    Parameters:
        ic_type(str):
            * '' - normal IC
            * 'T'
            * 'S'
            * 'C'
            * 'B'
            * 'M'
            * 'ZM'

        order(int): order(number of edge) of subgraph
    '''

    @classmethod
    def preset(cls):
        return (
            cls(t, o)
            for t in ['', 'T', 'S', 'C', 'B', 'M', 'ZM']
            for o in range(6)
        )

    @property
    def dependencies(self):
        if self.ic_type in ['T', 'S', 'C', 'B']:
            return dict(
                ICm=self.make_key('', self.order)
            )
        else:
            return dict(
                iAgs=Ag.make_key(self.order)
            )

    def __init__(self, ic_type='', order=0):
        assert ic_type in ['', 'T', 'S', 'C', 'B', 'M', 'ZM']
        self.ic_type = ic_type
        self.order = order

    @property
    def descriptor_name(self):
        return '{}IC{}'.format(self.ic_type, self.order)

    @property
    def descriptor_key(self):
        return self.make_key(self.ic_type, self.order)

    def calculate(self, mol, iAgs=None, ICm=None):
        N = mol.GetNumAtoms()
        if self.ic_type in ['', 'M', 'ZM']:
            ids = iAgs[0]
            Ags = iAgs[1]

            if self.ic_type == '':
                w = 1
            elif self.ic_type == 'M':
                w = np.vectorize(lambda i: mol.GetAtomWithIdx(int(i)).GetMass())(ids)
            else:
                w = Ags * np.vectorize(lambda i: mol.GetAtomWithIdx(int(i)).GetAtomicNum())(ids)

            e = np.vectorize(lambda Ag: entropy_term(float(Ag) / N))(Ags)
            return -np.sum(w * e)

        if self.ic_type == 'T':
            return N * ICm
        elif self.ic_type == 'S':
            return ICm / math.log(N, 2)
        elif self.ic_type == 'C':
            return math.log(N, 2) - ICm
        elif self.ic_type == 'B':
            bts = sum(b.GetBondTypeAsDouble() for b in mol.GetBonds())
            return ICm / math.log(bts, 2)


def entropy_term(a):
    return a * math.log(a, 2)
