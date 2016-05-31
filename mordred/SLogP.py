r"""Wildman-Crippen LogP/MR descriptor.

References
    * :cite:`10.1021/ci990307l`
"""
from rdkit.Chem import Crippen

from ._base import Descriptor


__all__ = ('SLogP', 'SMR',)


class WildmanCrippenBase(Descriptor):
    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return self.__class__.__name__

    def __reduce_ex__(self, version):
        return self.__class__, ()

    explicit_hydrogens = False

    rtype = float


class SLogP(WildmanCrippenBase):
    r"""Wildman-Crippen LogP descriptor(rdkit wrapper)."""

    __slots__ = ()

    def calculate(self, mol):
        return Crippen.MolLogP(mol)


class SMR(WildmanCrippenBase):
    r"""Wildman-Crippen MR descriptor(rdkit wrapper)."""

    __slots__ = ()

    def calculate(self, mol):
        return Crippen.MolMR(mol)


if __name__ == '__main__':
    from .__main__ import submodule
    submodule()
