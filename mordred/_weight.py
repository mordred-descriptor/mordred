from ._base import Descriptor
from rdkit.Chem.Descriptors import ExactMolWt


class Weight(Descriptor):
    r'''
    molecular weight descriptor

    Parameters:
        averaged(bool): averaged by number of atom

    Returns:
        float: exact molecular weight
    '''

    explicit_hydrogens = True
    require_connected = False

    @classmethod
    def preset(cls):
        yield cls(False)
        yield cls(True)

    def __str__(self):
        return 'AMW' if self.averaged else 'MW'

    descriptor_keys = 'averaged',

    def __init__(self, averaged=False):
        self.averaged = averaged

    def calculate(self, mol):
        w = ExactMolWt(mol)
        if self.averaged:
            w /= mol.GetNumAtoms()

        return w
