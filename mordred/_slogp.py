from rdkit.Chem import Crippen

from ._base import Descriptor


class WildmanCrippenBase(Descriptor):
    @classmethod
    def preset(cls):
        yield cls()

    def __str__(self):
        return self.__class__.__name__

    explicit_hydrogens = False


class SLogP(WildmanCrippenBase):
    r"""Wildman-Crippen LogP descriptor(rdkit wrapper).

    :rtype: float
    """

    __slots__ = ()

    def calculate(self, mol):
        return Crippen.MolLogP(mol)


class SMR(WildmanCrippenBase):
    r"""Wildman-Crippen MR descriptor(rdkit wrapper).

    :rtype: float
    """

    __slots__ = ()

    def calculate(self, mol):
        return Crippen.MolMR(mol)
