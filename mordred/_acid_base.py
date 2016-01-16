from ._base import Descriptor
from rdkit import Chem


class SmartsCount(Descriptor):
    def get_smarts(self):
        self.mols = [Chem.MolFromSmarts(s) for s in self.smarts]
        return self.mols

    def calculate(self, mol):
        mols = getattr(self, 'mols', None) or self.get_smarts()
        return sum(len(mol.GetSubstructMatches(q)) for q in mols)


class AcidicGroupCount(SmartsCount):
    r'''
    acidic group count descriptor

    Returns:
        int: number of acidic groups
    '''

    def __str__(self):
        return 'nAcid'

    smarts = [
        "[$([O;H1]-[C,S,P]=O)]",
        "[$([*;-;!$(*~[*;+])])]",
        "[$([NH](S(=O)=O)C(F)(F)F)]",
        "[$(n1nnnc1)]",
    ]


class BasicGroupCount(SmartsCount):
    r'''
    basic group count descriptor

    Returns:
        int: number of basic groups
    '''

    def __str__(self):
        return 'nBase'

    smarts = [
        "[$([NH2]-[CX4])]",
        "[$([NH](-[CX4])-[CX4])]",
        "[$(N(-[CX4])(-[CX4])-[CX4])]",
        "[$([*;+;!$(*~[*;-])])]",
        "[$(N=C-N)]",
        "[$(N-C=N)]",
    ]
