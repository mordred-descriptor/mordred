from ._base import Descriptor
from rdkit import Chem


class SmartsCount(Descriptor):
    def __init__(self):
        self.mols = []
        for s in self.smarts:
            mol = Chem.MolFromSmarts(s)
            if mol is None:
                raise ValueError('invalid smarts: {}'.format(s))
            self.mols.append(mol)

    def calculate(self, mol):
        return sum(len(mol.GetSubstructMatches(q)) for q in self.mols)


class AcidicGroupCount(SmartsCount):
    r'''
    acidic group count descriptor

    Returns:
        int: number of acidic groups
    '''

    descriptor_name = 'nAcid'
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

    descriptor_name = 'nBase'
    smarts = [
        "[$([NH2]-[CX4])]",
        "[$([NH](-[CX4])-[CX4])]",
        "[$(N(-[CX4])(-[CX4])-[CX4])]",
        "[$([*;+;!$(*~[*;-])])]",
        "[$(N=C-N)]",
        "[$(N-C=N)]",
    ]
