from .._base import Descriptor
from .._atomic_property import Polarizabilities78, Polarizabilities94

class APol(Descriptor):
    descriptor_name = 'apol'

    @property
    def descriptor_key(self):
        return self.make_key(self.use78)

    def __init__(self, use78=False):
        self.use78 = use78
        self.table = Polarizabilities78 if use78 else Polarizabilities94

    def calculate(self, mol):
        return sum((self.table[a.GetAtomicNum()] for a in mol.GetAtoms()))


class BPol(Descriptor):
    descriptor_name = 'bpol'

    @property
    def descriptor_key(self):
        return self.make_key(self.use78)

    def __init__(self, use78=False):
        self.use78 = use78
        self.table = Polarizabilities78 if use78 else Polarizabilities94

    def calculate(self, mol):
        def bond_pol(bond):
            a = bond.GetBeginAtom().GetAtomicNum()
            b = bond.GetEndAtom().GetAtomicNum()
            return abs(self.table[a] - self.table[b])

        return sum(bond_pol(b) for b in mol.GetBonds())


_descriptors = [APol, BPol]
__all__ = [d.__name__ for d in _descriptors]
