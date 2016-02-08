from rdkit import Chem


class MordredException(Exception):
    pass


class DescriptorException(MordredException):
    def __init__(self, desc, e, mol, parent=None):
        self.desc = desc
        self.e = e
        self.mol = mol
        self.parent = parent

    def __reduce_ex__(self, version):
        return self.__class__, (self.desc, self.e, self.mol, self.parent)

    def __str__(self):
        if self.parent is None:
            return '{}({!r}): {!r}'.format(
                self.desc,
                Chem.MolToSmiles(self.mol),
                self.e,
            )

        return '{}/{}({!r}): {!r}'.format(
            self.parent,
            self.desc,
            Chem.MolToSmiles(self.mol),
            self.e,
        )
