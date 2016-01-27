import math

from itertools import chain

from rdkit import Chem

from ._base import Descriptor


class PathCount(Descriptor):
    r"""path count descriptor.

    :type order: int
    :param order: path order(number of bonds in path)

    :type pi: bool
    :param pi: calculate pi-path count

    :type total: bool
    :param total: total path count(1 to order)

    :type log: bool
    :param log: use log scale

    :rtype: int(path-count) or float(pi-path-count)
    """

    explicit_hydrogens = False

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

    def __init__(self, order=1, pi=False, total=False, log=False):
        assert order >= 0

        self._order = order
        self._pi = pi
        self._total = total
        self._log = log

    def dependencies(self):
        if self._total:
            if self._order == 1:
                return dict(acc=self.__class__(
                    0, pi=self._pi, total=False
                ))

            else:
                return dict(acc=self.__class__(
                    self._order - 1,
                    self._pi,
                    self._total,
                ))

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

    def _find_paths(self, mol):
        for path in Chem.FindAllPathsOfLengthN(mol, self._order):
            aids = set()
            before = None
            w = 1

            for i in self._bond_ids_to_atom_ids(mol, path):
                if i in aids:
                    break

                aids.add(i)

                if self._pi and before is not None:
                    bond = mol.GetBondBetweenAtoms(before, i)
                    w *= bond.GetBondTypeAsDouble()

                before = i

            else:
                yield w

    def calculate(self, mol, acc=None):
        if self._order == 0:
            return mol.GetNumAtoms()

        v = sum(self._find_paths(mol))

        if acc is not None:
            v = acc + v

        if self._log:
            v = math.log(v + 1)

        return v
