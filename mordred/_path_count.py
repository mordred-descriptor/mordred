import math

from itertools import chain

from rdkit import Chem

from ._base import Descriptor


class PathCountBase(Descriptor):
    explicit_hydrogens = False


class PathCountCache(PathCountBase):
    __slots__ = ('_order',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._order,)

    def __init__(self, order):
        self._order = order

    @staticmethod
    def _bond_ids_to_atom_ids(mol, p):
        it = iter(p)

        try:
            b0 = mol.GetBondWithIdx(next(it))
        except StopIteration:
            return []

        try:
            b1 = mol.GetBondWithIdx(next(it))
        except StopIteration:
            return [b0.GetBeginAtomIdx(), b0.GetEndAtomIdx()]

        a0f, a0t = b0.GetBeginAtomIdx(), b0.GetEndAtomIdx()
        a1f, a1t = b1.GetBeginAtomIdx(), b1.GetEndAtomIdx()

        if a0f in [a1f, a1t]:
            path = [a0t, a0f]
            current = a1f if a0f == a1t else a1t
        else:
            path = [a0f, a0t]
            current = a1f if a0t == a1t else a1t

        for i in it:
            bn = mol.GetBondWithIdx(i)

            anf, ant = bn.GetBeginAtomIdx(), bn.GetEndAtomIdx()

            path.append(current)

            if anf == current:
                current = ant
            else:
                current = anf

        path.append(current)
        return path

    def calculate(self, mol):
        l = 0
        pi = 0

        for path in Chem.FindAllPathsOfLengthN(mol, self._order):
            aids = set()
            before = None
            w = 1

            for i in self._bond_ids_to_atom_ids(mol, path):
                if i in aids:
                    break

                aids.add(i)

                if before is not None:
                    bond = mol.GetBondBetweenAtoms(before, i)
                    w *= bond.GetBondTypeAsDouble()

                before = i

            else:
                l += 1
                pi += w

        return l, pi


class PathCount(PathCountBase):
    r"""path count descriptor.

    :type order: int
    :param order: path order(number of bonds in path)

    :type pi: bool
    :param pi: calculate pi-path count

    :type total: bool
    :param total: total path count(1 to order)

    :type log: bool
    :param log: use log scale
    """

    @classmethod
    def preset(cls):
        return chain(
            (cls(o, False, False, False) for o in range(2, 11)),
            [cls(10, False, True, False)],
            (cls(o, True, False, True) for o in range(1, 11)),
            [cls(10, True, True, True)],
        )

    def __str__(self):
        if self._total:
            base = 'T'
        else:
            base = ''

        if self._pi:
            base += 'piPC'
        else:
            base += 'MPC'

        return base + str(self._order)

    __slots__ = ('_order', '_pi', '_total', '_log',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._order, self._pi, self._total, self._log)

    def __init__(self, order=1, pi=False, total=False, log=False):
        assert order >= 0

        self._order = order
        self._pi = pi
        self._total = total
        self._log = log

    def dependencies(self):
        deps = {'PC': PathCountCache(self._order)}
        if self._total and self._order > 0:
            deps['acc'] = self.__class__(
                self._order - 1,
                self._pi,
                self._total,
            )

        return deps

    def calculate(self, mol, PC, acc=None):
        if self._order == 0:
            return self.rtype(mol.GetNumAtoms())

        v = PC[1] if self._pi else PC[0]

        if acc is not None:
            v = acc + v

        if self._log:
            v = math.log(v + 1)

        return v

    @property
    def rtype(self):
        r"""
        * pi = True: :py:class:`float`
        * pi = False: :py:class:`int`
        """
        return float if self._pi else int
