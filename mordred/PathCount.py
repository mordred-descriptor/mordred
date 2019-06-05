import math
from itertools import chain

from rdkit import Chem

from ._base import Descriptor

__all__ = ("PathCount",)


class PathCountBase(Descriptor):
    __slots__ = ()
    explicit_hydrogens = False


class PathCountCache(PathCountBase):
    __slots__ = ("_order", "_bonds")

    def parameters(self):
        return (self._order,)

    def __init__(self, order):
        self._order = order

    def _gen_bonds(self):
        self._bonds = [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in self.mol.GetBonds()
        ]

    def _bond_ids_to_atom_ids(self, p):
        it = iter(p)

        try:
            a0f, a0t = self._bonds[next(it)]
        except StopIteration:
            return []

        try:
            a1f, a1t = self._bonds[next(it)]
        except StopIteration:
            return a0f, a0t

        if a0f in [a1f, a1t]:
            path = [a0t, a0f]
            current = a1f if a0f == a1t else a1t
        else:
            path = [a0f, a0t]
            current = a1f if a0t == a1t else a1t

        for i in it:
            anf, ant = self._bonds[i]

            path.append(current)

            if anf == current:
                current = ant
            else:
                current = anf

        path.append(current)
        return path

    def calculate(self):
        L = 0
        pi = 0

        self._gen_bonds()

        for path in Chem.FindAllPathsOfLengthN(self.mol, self._order):
            aids = set()
            before = None
            w = 1

            for i in self._bond_ids_to_atom_ids(path):
                if i in aids:
                    break

                aids.add(i)

                if before is not None:
                    bond = self.mol.GetBondBetweenAtoms(before, i)
                    w *= bond.GetBondTypeAsDouble()

                before = i

            else:
                L += 1
                pi += w

        return L, pi


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

    since = "1.0.0"
    __slots__ = ("_order", "_pi", "_total", "_log")

    def description(self):
        return "{}-ordered {}{}path count{}".format(
            self._order,
            "total " if self._total else "",
            "pi-" if self._pi else "",
            " (log scale)" if self._log else "",
        )

    @classmethod
    def preset(cls, version):
        return chain(
            (cls(o, False, False, False) for o in range(2, 11)),
            [cls(10, False, True, False)],
            (cls(o, True, False, True) for o in range(1, 11)),
            [cls(10, True, True, True)],
        )

    def __str__(self):
        base = "T" if self._total else ""
        pi = "piPC" if self._pi else "MPC"

        return "{}{}{}".format(base, pi, self._order)

    def parameters(self):
        return self._order, self._pi, self._total, self._log

    def __init__(self, order=1, pi=False, total=False, log=False):
        assert order >= 0

        self._order = order
        self._pi = pi
        self._total = total
        self._log = log

    def dependencies(self):
        deps = {"PC": PathCountCache(self._order)}
        if self._total and self._order > 0:
            deps["acc"] = self.__class__(self._order - 1, self._pi, self._total)

        return deps

    def calculate(self, PC, acc=None):
        if self._order == 0:
            return self.rtype(self.mol.GetNumAtoms())

        v = PC[1] if self._pi else PC[0]

        if acc is not None:
            v = acc + v

        if self._log:
            v = math.log(v + 1)

        return v

    @property
    def rtype(self):
        r"""Return type.

        * pi = True: :py:class:`float`
        * pi = False: :py:class:`int`
        """
        return float if self._pi else int
