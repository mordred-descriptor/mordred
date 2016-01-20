from rdkit.Chem import Crippen as _Crippen

from ._base import Descriptor


class WildmanCrippenLogP(Descriptor):
    r"""Wildman-Crippen LogP/MR descriptor.

    :type prop: str
    :param type: 'LogP' or 'MR'

    :rtype: float

    References
        * :cite:`10.1021/ci990307l`
    """

    @classmethod
    def preset(cls):
        yield cls('LogP')
        yield cls('MR')

    explicit_hydrogens = False
    require_connected = False

    def __str__(self):
        return 'Crippen{}'.format(self.prop)

    __slots__ = ('prop',)

    def __init__(self, prop='LogP'):
        assert prop in ['LogP', 'MR']
        self.prop = prop

    def calculate(self, mol):
        if self.prop == 'LogP':
            return _Crippen.MolLogP(mol)
        else:
            return _Crippen.MolMR(mol)
