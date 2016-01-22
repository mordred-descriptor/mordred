from abc import abstractproperty

from rdkit import Chem

from ._base import Descriptor


class SmartsCountBase(Descriptor):
    require_connected = False

    @classmethod
    def preset(cls):
        yield cls()

    def _get_smarts(self):
        self._mols = [Chem.MolFromSmarts(s) for s in self.SMARTS]
        return self._mols

    def calculate(self, mol):
        mols = getattr(self, 'mols', None) or self._get_smarts()
        return sum(len(mol.GetSubstructMatches(q)) for q in mols)

    @abstractproperty
    def SMARTS(self):
        u"""target smarts.

        (abstruct property)
        """
        pass


class AcidicGroupCount(SmartsCountBase):
    r"""acidic group count descriptor.

    :rtype: int
    """

    __slots__ = ()

    def __str__(self):
        return 'nAcid'

    SMARTS = (
        "[$([O;H1]-[C,S,P]=O)]",
        "[$([*;-;!$(*~[*;+])])]",
        "[$([NH](S(=O)=O)C(F)(F)F)]",
        "[$(n1nnnc1)]",
    )


class BasicGroupCount(SmartsCountBase):
    r"""basic group count descriptor.

    :rtype: int
    """

    __slots__ = ()

    def __str__(self):
        return 'nBase'

    SMARTS = (
        "[$([NH2]-[CX4])]",
        "[$([NH](-[CX4])-[CX4])]",
        "[$(N(-[CX4])(-[CX4])-[CX4])]",
        "[$([*;+;!$(*~[*;-])])]",
        "[$(N=C-N)]",
        "[$(N-C=N)]",
    )
