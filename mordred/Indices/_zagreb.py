from .._base import Descriptor
from .._common import Valence


class ZagrebIndex(Descriptor):
    explicit_hydrogens = False

    descriptor_defaults = [(1, 1), (2, 1), (1, -1), (2, -1)]

    @property
    def dependencies(self):
        return dict(
            V=Valence.make_key(
                self.explicit_hydrogens,
                False,
            ),
        )

    @property
    def descriptor_name(self):
        if self.variable in [1, -1]:
            m = '' if self.variable == 1 else 'm'
            return '{}Zagreb{}'.format(m, self.version)

        return 'Zagreb{}_{}'.format(self.version, self.variable)

    def __init__(self, version=1, variable=1):
        assert version in [1, 2]
        self.version = version
        self.variable = variable

    @property
    def descriptor_key(self):
        return self.make_key(self.version, self.variable)

    def calculate(self, mol, V):
        if not isinstance(self.variable, int) or self.variable < 0:
            V = V.astype('float')

        if self.version == 1:
            return (V ** (self.variable * 2)).sum()
        else:
            return sum(
                (V[b.GetBeginAtomIdx()] * V[b.GetEndAtomIdx()]) ** self.variable
                for b in mol.GetBonds()
            )
