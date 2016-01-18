from rdkit.Chem.Descriptors import ExactMolWt

from ._base import Descriptor


class Weight(Descriptor):
    r"""molecular weight descriptor.

    :type averaged: bool
    :param averaged: averaged by number of atom

    :rtype: float
    """

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
