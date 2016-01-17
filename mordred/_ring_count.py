from ._base import Descriptor
from rdkit import Chem
import networkx


class RingCountBase(Descriptor):
    explicit_hydrogens = False
    require_connected = False


class Rings(RingCountBase):
    def calculate(self, mol):
        return [frozenset(s) for s in Chem.GetSymmSSSR(mol)]


class FusedRings(RingCountBase):
    @property
    def dependencies(self):
        return dict(Rings=Rings())

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
        order(int or None): number of bonds in ring
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

    def __str__(self):
        attrs = []

        if self.greater:
            attrs.append('G')

        if self.order is not None:
            attrs.append(str(self.order))

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

    descriptor_keys = 'order', 'greater', 'fused', 'aromatic', 'hetero'

    def __init__(self, order=None, greater=False, fused=False, aromatic=None, hetero=None):
        self.order = order
        self.greater = greater
        self.fused = fused
        self.aromatic = aromatic
        self.hetero = hetero

    @property
    def dependencies(self):
        return dict(
            Rs=(FusedRings if self.fused else Rings)()
        )

    def _check_order(self, R):
        if self.order is None:
            return True

        if self.greater:
            return len(R) >= self.order
        else:
            return len(R) == self.order

    def _check_arom(self, mol, R):
        if self.aromatic is None:
            return True

        is_arom = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in R)

        if self.aromatic:
            return is_arom

        return not is_arom

    def _check_hetero(self, mol, R):
        if self.hetero is None:
            return True

        has_hetero = any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in R)

        if self.hetero:
            return has_hetero

        return not has_hetero

    def calculate(self, mol, Rs):
        return sum(
            1 for R in Rs
            if self._check_order(R) and self._check_arom(mol, R) and self._check_hetero(mol, R)
        )
