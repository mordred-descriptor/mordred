import numpy as np

from rdkit import Chem

from scipy.sparse.csgraph import shortest_path

from . import _atomic_property
from ._base import Descriptor
from ._matrix_attributes import get_method, methods


class BaryszMatrixBase(Descriptor):
    explicit_hydrogens = False
    __slots__ = ()

    @property
    def gasteiger_charges(self):
        return getattr(self._prop, 'gasteiger_charges', False)


_carbon = Chem.Atom(6)


class Barysz(BaryszMatrixBase):
    __slots__ = ('_prop',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._prop,)

    def __init__(self, prop):
        self._prop = prop

    def calculate(self, mol):
        N = mol.GetNumAtoms()

        C = self._prop(_carbon)

        dmat = np.zeros((N, N))

        for bond in mol.GetBonds():
            ai = bond.GetBeginAtom()
            aj = bond.GetEndAtom()

            i = ai.GetIdx()
            j = aj.GetIdx()

            pi = bond.GetBondTypeAsDouble()

            w = float(C * C) / float(self._prop(ai) * self._prop(aj) * pi)
            dmat[i, j] = w
            dmat[j, i] = w

        if np.any(~np.isfinite(dmat)):
            return None

        sp = shortest_path(dmat)
        np.fill_diagonal(sp, [1. - float(C) / self._prop(a) for a in mol.GetAtoms()])
        return sp


class BaryszMatrix(BaryszMatrixBase):
    r"""barysz matrix descriptor.

    :type prop: :py:class:`str` or :py:class:`function`
    :param prop: :ref:`atomic_properties`

    :type type: str
    :param type: :ref:`matrix_aggregating_methods`

    :returns: NaN when any properties are NaN
    """

    @classmethod
    def preset(cls):
        return (cls(p, m) for p in _atomic_property.get_properties() for m in methods)

    def __str__(self):
        return '{}_Dz{}'.format(self._type.__name__, self._prop_name)

    __slots__ = ('_prop_name', '_prop', '_type',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._prop, self._type)

    def __init__(self, prop='Z', type='SpMax'):
        self._prop_name, self._prop = _atomic_property.getter(prop, self.explicit_hydrogens)
        self._type = get_method(type)

    def dependencies(self):
        return dict(
            result=self._type(
                Barysz(self._prop),
                self.explicit_hydrogens,
                self.gasteiger_charges,
                self.kekulize,
            )
        )

    def calculate(self, mol, result):
        return result

    rtype = float
