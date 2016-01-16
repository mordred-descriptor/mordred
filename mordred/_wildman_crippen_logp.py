from ._base import Descriptor
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

    def __str__(self):
        return 'Crippen{}'.format(self.value)

    descriptor_keys = 'value',

    def __init__(self, value='LogP'):
        assert value in ['LogP', 'MR']
        self.value = value

    def calculate(self, mol):
        if self.value == 'LogP':
            return _Crippen.MolLogP(mol)
        else:
            return _Crippen.MolMR(mol)
