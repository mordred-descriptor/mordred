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

    def __str__(self):
        return self._name + str(self._order)

    @classmethod
    def preset(cls):
        return (cls(o) for o in range(6))

    __slots__ = ('_order',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._order,)

    def __init__(self, order=0):
        self._order = order

    rtype = float


class Ag(InformationContentBase):
    @classmethod
    def preset(cls):
        return ()

    __slots__ = ('_order',)

    def dependencies(self):
        return dict(
            D=DistanceMatrix(self.explicit_hydrogens)
        )

    def calculate(self, mol, D):
        atoms = [neighborhood_code(mol, i, self._order) for i in range(mol.GetNumAtoms())]
        ad = {a: i for i, a in enumerate(atoms)}
        Ags = [(k, sum(1 for _ in g)) for k, g in groupby(sorted(atoms))]
        return (
            np.array([ad[k] for k, _ in Ags]),
            np.array([ag for _, ag in Ags], dtype='float'),
        )


def _shannon_entropy_term(a):
    return a * np.log2(a)


shannon_entropy_term = np.vectorize(_shannon_entropy_term)


def shannon_entropy(a, w=1):
    N = np.sum(a)
    return -np.sum(w * shannon_entropy_term(a / N))


class InformationContent(InformationContentBase):
    r"""neighborhood information content descriptor.

    :type order: int
    :param order: order(number of edge) of subgraph
    """

    _name = 'IC'

    def dependencies(self):
        return dict(iAgs=Ag(self._order))

    def calculate(self, mol, iAgs):
        _, Ags = iAgs
        return shannon_entropy(Ags)


class TotalIC(InformationContentBase):
    r"""neighborhood total information content descriptor.

    .. math::
        {\rm TIC}_m = A \cdot {\rm IC}_m

    :type order: int
    :param order: order(number of edge) of subgraph
    """

    _name = 'TIC'

    def dependencies(self):
        return dict(ICm=InformationContent(self._order))

    def calculate(self, mol, ICm):
        A = mol.GetNumAtoms()

        return A * ICm


class StructuralIC(InformationContentBase):
    r"""structural information content descriptor.

    .. math::
        {\rm SIC}_m = \frac{{\rm IC}_m}{\log_2 A}

    :type order: int
    :param order: order(number of edge) of subgraph
    """

    _name = 'SIC'

    def dependencies(self):
        return dict(ICm=InformationContent(self._order))

    def calculate(self, mol, ICm):
        A = mol.GetNumAtoms()

        return ICm / np.log2(A)


class BondingIC(InformationContentBase):
    r"""bonding information content descriptor.

    .. math::
        {\rm BIC}_m = \frac{{\rm IC}_m}{\log_2 \sum^B_{b=1} \pi^{*}_b}

    :type order: int
    :param order: order(number of edge) of subgraph

    :returns: NaN when :math:`\sum^B_{b=1} \pi^{*}_b <= 0`
    """

    _name = 'BIC'

    def dependencies(self):
        return dict(ICm=InformationContent(self._order))

    def calculate(self, mol, ICm):
        B = sum(b.GetBondTypeAsDouble() for b in mol.GetBonds())

        if B == 0:
            return np.nan

        log2B = np.log2(B)
        if log2B == 0:
            return np.nan

        return ICm / log2B


class ComplementaryIC(InformationContentBase):
    r"""complementary information content descriptor.

    .. math::
        {\rm CIC}_m = \log_2 A - {\rm IC}_m

    :type order: int
    :param order: order(number of edge) of subgraph
    """

    _name = 'CIC'

    def dependencies(self):
        return dict(ICm=InformationContent(self._order))

    def calculate(self, mol, ICm):
        A = mol.GetNumAtoms()

        return np.log2(A) - ICm


class ModifiedIC(InformationContentBase):
    r"""modified information content index descriptor.

    :type order: int
    :param order: order(number of edge) of subgraph
    """

    _name = 'MIC'

    def dependencies(self):
        return dict(iAgs=Ag(self._order))

    def calculate(self, mol, iAgs):
        ids, Ags = iAgs
        w = np.vectorize(lambda i: mol.GetAtomWithIdx(int(i)).GetMass())(ids)
        return shannon_entropy(Ags, w)


class ZModifiedIC(InformationContentBase):
    r"""Z-modified information content index descriptor.

    :type order: int
    :param order: order(number of edge) of subgraph
    """

    _name = 'ZMIC'

    def dependencies(self):
        return dict(iAgs=Ag(self._order))

    def calculate(self, mol, iAgs):
        ids, Ags = iAgs
        w = Ags * np.vectorize(lambda i: mol.GetAtomWithIdx(int(i)).GetAtomicNum())(ids)
        return shannon_entropy(Ags, w)
