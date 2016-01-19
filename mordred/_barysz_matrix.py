import numpy as np

from rdkit import Chem

from scipy.sparse.csgraph import shortest_path

from . import _atomic_property
from ._base import Descriptor
from ._matrix_attributes import get_method, methods


class BaryszMatrixBase(Descriptor):
    explicit_hydrogens = False

    @property
    def gasteiger_charges(self):
        return getattr(self.prop, 'gasteiger_charges', False)


_carbon = Chem.Atom(6)


class Barysz(BaryszMatrixBase):
    descriptor_keys = 'prop',

    def __init__(self, prop):
        self.prop = prop

    def calculate(self, mol):
        N = mol.GetNumAtoms()

        C = self.prop(_carbon)

        dmat = np.zeros((N, N))

        for bond in mol.GetBonds():
            ai = bond.GetBeginAtom()
            aj = bond.GetEndAtom()

            i = ai.GetIdx()
            j = aj.GetIdx()

            pi = bond.GetBondTypeAsDouble()

            w = float(C * C) / float(self.prop(ai) * self.prop(aj) * pi)
            dmat[i, j] = w
            dmat[j, i] = w

        sp = shortest_path(dmat)
        np.fill_diagonal(sp, [1. - float(C) / self.prop(a) for a in mol.GetAtoms()])
        return sp


class BaryszMatrix(BaryszMatrixBase):
    r"""barysz matrix descriptor.

    :type prop: str or function
    :param prop: :ref:`atomic_properties`

    :type type: str
    :param type: :ref:`matrix_aggregating_methods`

    :rtype: float
    """

    @classmethod
    def preset(cls):
        return (cls(p, m) for p in _atomic_property.get_properties() for m in methods)

    def __str__(self):
        return '{}_Dz{}'.format(self.type.__name__, self.prop_name)

    descriptor_keys = 'prop', 'type'

    def __init__(self, prop='Z', type='SpMax'):
        self.prop_name, self.prop = _atomic_property.getter(prop, self.explicit_hydrogens)
        self.type = get_method(type)

    def dependencies(self):
        return dict(
            result=self.type(
                Barysz(self.prop),
                self.explicit_hydrogens,
                self.gasteiger_charges,
                self.kekulize,
            )
        )

    def calculate(self, mol, result):
        return result
