from ._base import Descriptor


class AtomCount(Descriptor):
    r'''
    atom count descriptor

    Parameters:
        type(str): type to count. 'Atom', 'HeavyAtom', 'X'(all halogen), or element symbol.

    Returns:
        int: count
    '''

    @classmethod
    def preset(cls):
        return map(cls, [
            'Atom', 'HeavyAtom',
            'H', 'B', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'X',
        ])

    @property
    def explicit_hydrogens(self):
        return self.type in set(['H', 'Atom'])

    def __str__(self):
        return 'n' + self.type

    descriptor_keys = 'type',

    def __init__(self, type='Atom'):
        self.type = type

    def _calc_X(self, mol):
        X = set([9, 17, 35, 53, 85, 117])
        return sum(a.GetAtomicNum() in X for a in mol.GetAtoms())

    def _calc(self, mol):
        return sum(a.GetSymbol() == self.type for a in mol.GetAtoms())

    def _calc_all(self, mol):
        return mol.GetNumAtoms()

    def calculate(self, mol):
        if self.type == 'X':
            return self._calc_X(mol)
        elif self.type in ['Atom', 'HeavyAtom']:
            return self._calc_all(mol)
        else:
            return self._calc(mol)
