from .._base import Descriptor
from rdkit.Chem import Crippen as _Crippen


class WildmanCrippenLogP(Descriptor):
    r'''
    Wildman-Crippen LogP/MR descriptor

    Parameters:
        value(str): 'LogP' or 'MR'

    Returns:
        float: LogP or MR value
    '''

    @classmethod
    def preset(cls):
        yield cls('LogP')
        yield cls('MR')

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

_descriptors = [WildmanCrippenLogP]
__all__ = [d.__name__ for d in _descriptors]
