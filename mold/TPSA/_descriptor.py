from .._base import Descriptor
from rdkit.Chem import rdMolDescriptors


class TPSA(Descriptor):
    descriptor_defaults = [()]

    @property
    def descriptor_name(self):
        return 'TPSA'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol):
        return rdMolDescriptors.CalcTPSA(mol)

_descriptors = [TPSA]
__all__ = [d.__name__ for d in _descriptors]
