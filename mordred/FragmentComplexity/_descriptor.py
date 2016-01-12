from .._base import Descriptor


class FragmentComplexity(Descriptor):
    explicit_hydrogens = False

    @property
    def descriptor_name(self):
        return 'fragCpx'

    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol):
        A = mol.GetNumAtoms()
        B = mol.GetNumBonds()
        H = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 6)
        return abs(B ** 2 - A ** 2 + A) + float(H) / 100


_descriptors = [FragmentComplexity]
__all__ = [d.__name__ for d in _descriptors]
