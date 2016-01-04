from .._base import Descriptor
from rdkit.Chem import Crippen as _Crippen


class WildmanCrippen(Descriptor):
    descriptor_defaults = [('LogP',), ('MR',)]

    explicit_hydrogens = False

    def __init__(self, value='LogP'):
        assert value in ['LogP', 'MR']
        self.value = value

    @property
    def descriptor_name(self):
        return 'Crippen{}'.format(self.value)

    @property
    def descriptor_key(self):
        return self.make_key(self.value)

    def calculate(self, mol):
        if self.value == 'LogP':
            return _Crippen.MolLogP(mol)
        else:
            return _Crippen.MolMR(mol)

_descriptors = [WildmanCrippen]
__all__ = [d.__name__ for d in _descriptors]
