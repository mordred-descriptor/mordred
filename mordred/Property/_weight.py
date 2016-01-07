from .._base import Descriptor
from rdkit.Chem.Descriptors import ExactMolWt


class Weight(Descriptor):
    explicit_hydrogens = True

    descriptor_defaults = [(False,), (True,)]

    @property
    def descriptor_name(self):
        return 'AMW' if self.averaged else 'MW'

    def __init__(self, averaged=False):
        self.averaged = averaged

    @property
    def descriptor_key(self):
        return self.make_key(self.averaged)

    def calculate(self, mol):
        w = ExactMolWt(mol)
        if self.averaged:
            w /= mol.GetNumAtoms()

        return w
