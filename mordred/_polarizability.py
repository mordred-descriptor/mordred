from ._base import Descriptor
from ._atomic_property import Polarizabilities78, Polarizabilities94


class APol(Descriptor):
    r'''
    atomic polarizability descriptor

    Parameters:
        use78(bool): use old atomic polarizability data

    Returns:
        float: sum of atomic polarizability
    '''

    def __str__(self):
        return 'apol78' if self.use78 else 'apol'

    descriptor_keys = 'use78',

    def __init__(self, use78=False):
        self.use78 = use78
        self.table = Polarizabilities78 if use78 else Polarizabilities94

    def calculate(self, mol):
        return sum(self.table[a.GetAtomicNum()] for a in mol.GetAtoms())


class BPol(Descriptor):
    r'''
    bond polarizability descriptor

    Parameters:
        use78(bool): use old atomic polarizability data

    Returns:
        float: sum of bond polarizability
    '''

    def __str__(self):
        return 'bpol78' if self.use78 else 'bpol'

    descriptor_keys = 'use78',

    def __init__(self, use78=False):
        self.use78 = use78
        self.table = Polarizabilities78 if use78 else Polarizabilities94

    def calculate(self, mol):
        def bond_pol(bond):
            a = bond.GetBeginAtom().GetAtomicNum()
            b = bond.GetEndAtom().GetAtomicNum()
            return abs(self.table[a] - self.table[b])

        return sum(bond_pol(b) for b in mol.GetBonds())
