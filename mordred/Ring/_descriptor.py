from .._base import Descriptor
from rdkit import Chem
import networkx


class RingCountBase(Descriptor):
    explicit_hydrogens = False


class Rings(RingCountBase):
    def calculate(self, mol):
        return [frozenset(s) for s in Chem.GetSymmSSSR(mol)]


class FusedRings(RingCountBase):
    @property
    def dependencies(self):
        return dict(Rings=Rings.make_key())

    def calculate(self, mol, Rings):
        if len(Rings) < 2:
            return []

        G = networkx.Graph()

        l = len(Rings)
        for i, j in ((i, j) for i in range(l) for j in range(i + 1, l)):
            if len(Rings[i] & Rings[j]) >= 2:
                G.add_edge(i, j)

        return [
            frozenset(j for i in ring_ids for j in Rings[i])
            for ring_ids in networkx.connected_components(G)
        ]


class RingCount(RingCountBase):
    r'''
    ring count descriptor

    Parameters:
        length(int or None): number of bonds in ring
        greater(bool): count length or greater rings
        fused(bool): count fused rings
        aromatic(bool, None):
            * True: count aromatic rings
            * False: count non-aromatic rings
            * None: count any rings
        hetero(boo, None):
            * True: count hetero rings
            * False: count carbon rings
            * None: count any rings

    Returns:
        int: ring count
    '''

    @classmethod
    def preset(cls):
        for fused in [False, True]:
            for arom in [None, True, False]:
                for hetero in [None, True]:
                    yield cls(None, False, fused, arom, hetero)
                    for n in range(4 if fused else 3, 13):
                        yield cls(n, False, fused, arom, hetero)

                    yield cls(12, True, fused, arom, hetero)

    @property
    def descriptor_name(self):
        attrs = []

        if self.greater:
            attrs.append('G')

        if self.length is not None:
            attrs.append(str(self.length))

        if self.fused:
            attrs.append('F')

        if self.aromatic is True:
            attrs.append('a')
        elif self.aromatic is False:
            attrs.append('A')

        if self.hetero is True:
            attrs.append('H')
        elif self.hetero is False:
            attrs.append('C')

        return 'n{}Ring'.format(''.join(attrs))

    def __init__(self, length=None, greater=False, fused=False, aromatic=None, hetero=None):
        self.length = length
        self.greater = greater
        self.fused = fused
        self.aromatic = aromatic
        self.hetero = hetero

    @property
    def descriptor_key(self):
        return self.make_key(
            self.length,
            self.greater,
            self.fused,
            self.aromatic,
            self.hetero,
        )

    @property
    def dependencies(self):
        return dict(
            Rs=(FusedRings if self.fused else Rings).make_key()
        )

    def check_length(self, R):
        if self.length is None:
            return True

        if self.greater:
            return len(R) >= self.length
        else:
            return len(R) == self.length

    def check_arom(self, mol, R):
        if self.aromatic is None:
            return True

        is_arom = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in R)

        if self.aromatic:
            return is_arom

        return not is_arom

    def check_hetero(self, mol, R):
        if self.hetero is None:
            return True

        has_hetero = any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in R)

        if self.hetero:
            return has_hetero

        return not has_hetero

    def calculate(self, mol, Rs):
        return sum(
            1 for R in Rs
            if self.check_length(R) and self.check_arom(mol, R) and self.check_hetero(mol, R)
        )


_descriptors = [RingCount]
__all__ = [d.__name__ for d in _descriptors]
