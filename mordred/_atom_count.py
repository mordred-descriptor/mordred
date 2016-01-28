from rdkit.Chem import rdMolDescriptors

from ._base import Descriptor


class AtomCount(Descriptor):
    r"""atom count descriptor.

    :type type: str
    :param type: type to count.

        * 'Atom'
        * 'HeavyAtom'
        * 'Spiro'
        * 'Bridgehead'
        * 'X' - all halogen
        * element symbol
    """

    @classmethod
    def preset(cls):
        return map(cls, [
            'Atom', 'HeavyAtom', 'Spiro', 'Bridgehead',
            'H', 'B', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'X',
        ])

    @property
    def explicit_hydrogens(self):
        u"""require explicit_hydrogens when type is 'H' or 'Atom'."""
        return self._type in set(['H', 'Atom'])

    def __str__(self):
        return 'n' + self._type

    __slots__ = ('_type',)

    def __reduce_ex__(self, version):
        return self.__class__, (self._type,)

    def __init__(self, type='Atom'):
        self._type = type

    def _calc_X(self, mol):
        X = set([9, 17, 35, 53, 85, 117])
        return sum(a.GetAtomicNum() in X for a in mol.GetAtoms())

    def _calc(self, mol):
        return sum(a.GetSymbol() == self._type for a in mol.GetAtoms())

    def _calc_all(self, mol):
        return mol.GetNumAtoms()

    def calculate(self, mol):
        if self._type == 'X':
            return self._calc_X(mol)
        elif self._type in ['Atom', 'HeavyAtom']:
            return self._calc_all(mol)
        elif self._type == 'Spiro':
            return rdMolDescriptors.CalcNumSpiroAtoms(mol)
        elif self._type == 'Bridgehead':
            return rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        else:
            return self._calc(mol)

    rtype = int
