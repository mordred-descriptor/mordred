from ._base import Descriptor
from rdkit.Chem import Crippen as _Crippen


class WildmanCrippenLogP(Descriptor):
    r'''
    Wildman-Crippen LogP/MR descriptor

    Parameters:
        prop(str): 'LogP' or 'MR'

    Returns:
        float: LogP or MR value
    '''

    @classmethod
    def preset(cls):
        yield cls('LogP')
        yield cls('MR')

    explicit_hydrogens = False

    def __str__(self):
        return 'Crippen{}'.format(self.prop)

    descriptor_keys = 'prop',

    def __init__(self, prop='LogP'):
        assert prop in ['LogP', 'MR']
        self.prop = prop

    def calculate(self, mol):
        if self.prop == 'LogP':
            return _Crippen.MolLogP(mol)
        else:
            return _Crippen.MolMR(mol)
