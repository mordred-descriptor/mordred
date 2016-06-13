from abc import abstractproperty

from rdkit import Chem

from ._base import Descriptor

__all__ = ('AcidicGroupCount', 'BasicGroupCount',)


class SmartsCountBase(Descriptor):
    __slots__ = '_mol',

    @classmethod
    def preset(cls):
        yield cls()

    def _create_smarts(self):
        s = ','.join('$(' + s + ')' for s in self.SMARTS)
        self._mol = Chem.MolFromSmarts('[' + s + ']')
        return self._mol

    @abstractproperty
    def SMARTS(self):
        pass

    def __str__(self):
        return self._name

    def as_key(self):
        return self.__class__, ()

    def calculate(self):
        pat = getattr(self, '_mol', None) or self._create_smarts()
        return len(self.mol.GetSubstructMatches(pat))

    rtype = int


class AcidicGroupCount(SmartsCountBase):
    r"""acidic group count descriptor."""

    __slots__ = ()

    _name = 'nAcid'

    SMARTS = (
        '[O;H1]-[C,S,P]=O',
        '[*;-;!$(*~[*;+])]',
        '[NH](S(=O)=O)C(F)(F)F',
        'n1nnnc1',
    )


class BasicGroupCount(SmartsCountBase):
    r"""basic group count descriptor."""

    __slots__ = ()

    _name = 'nBase'

    SMARTS = (
        '[NH2]-[CX4]',
        '[NH](-[CX4])-[CX4]',
        'N(-[CX4])(-[CX4])-[CX4]',
        '[*;+;!$(*~[*;-])]',
        'N=C-N',
        'N-C=N',
    )
