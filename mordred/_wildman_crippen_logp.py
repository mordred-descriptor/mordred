from rdkit.Chem import Crippen as _Crippen

from ._base import Descriptor


class WildmanCrippenBase(Descriptor):
    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return self.__class__.__name__[7:]

    explicit_hydrogens = False
    require_connected = False


class WildmanCrippenLogP(WildmanCrippenBase):
    r"""Wildman-Crippen LogP descriptor.

    :rtype: float
    """

    __slots__ = ()

    def calculate(self, mol):
        return _Crippen.MolLogP(mol)


class WildmanCrippenMR(WildmanCrippenBase):
    r"""Wildman-Crippen MR descriptor.

    :rtype: float
    """

    __slots__ = ()

    def calculate(self, mol):
        return _Crippen.MolMR(mol)
