import numpy as np

from . import _atomic_property
from ._base import Descriptor


class BCUTBase(Descriptor):
    explicit_hydrogens = False
    require_connected = True


class Burden(BCUTBase):
    __slots__ = ()

    def __reduce_ex__(self, version):
        return self.__class__, ()

    def calculate(self, mol):
        N = mol.GetNumAtoms()

        mat = 0.001 * np.ones((N, N))

        for bond in mol.GetBonds():
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()
            i = a.GetIdx()
            j = b.GetIdx()

            try:
                w = bond.GetBondTypeAsDouble() / 10.0
            except RuntimeError:
                w = 1.0

            if a.GetDegree() == 1 or b.GetDegree() == 1:
                w += 0.01

            mat[i, j] = w
            mat[j, i] = w

        return mat


class BurdenEigenValues(BCUTBase):
    __slots__ = ('_prop', 'gasteiger_charges',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._prop, self.gasteiger_charges)

    def __init__(self, prop, gasteiger_charges):
        self._prop = prop
        self.gasteiger_charges = gasteiger_charges

    def dependencies(self):
        return dict(burden=Burden())

    def calculate(self, mol, burden):
        bmat = burden.copy()
        ps = np.array([self._prop(a) for a in mol.GetAtoms()])
        if np.any(np.isnan(ps)):
            return np.array([np.nan])

        np.fill_diagonal(bmat, ps)
        ev = np.linalg.eig(bmat)[0]

        if np.iscomplexobj(ev):
            ev = ev.real

        return np.sort(ev)[-1::-1]


class BCUT(BCUTBase):
    r"""BCUT descriptor.

    :type prop: :py:class:`str` or :py:class:`function`
    :param prop: :ref:`atomic_properties`

    :type nth: int
    :param nth: n-th eigen value. 0 is highest, -1 is lowest.

    :rtype: float
    :returns: NaN when

        * any atomic properties are NaN
        * :math:`\left| nth \right| > A`
    """

    @classmethod
    def preset(cls):
        return (
            cls(a, n)
            for a in _atomic_property.get_properties(istate=True, charge=True)
            for n in [0, -1]
        )

    @property
    def gasteiger_charges(self):
        return getattr(self._prop, 'gasteiger_charges', False)

    def __str__(self):
        if self._nth < 0:
            return 'BCUT{}-{}l'.format(self._prop_name, np.abs(self._nth))
        else:
            return 'BCUT{}-{}h'.format(self._prop_name, self._nth + 1)

    __slots__ = ('_prop', '_prop_name', '_nth',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._prop, self._nth)

    def __init__(self, prop='m', nth=0):
        self._prop_name, self._prop = _atomic_property.getter(prop, self.explicit_hydrogens)
        self._nth = nth

    def dependencies(self):
        return dict(bev=BurdenEigenValues(self._prop, self.gasteiger_charges))

    def calculate(self, mol, bev):
        try:
            return bev[self._nth]
        except IndexError:
            return np.nan
