from collections import defaultdict
from .._base import Descriptor


class CarbonTypesBase(Descriptor):
    explicit_hydrogens = False
    kekulize = True


class CarbonTypesCache(CarbonTypesBase):
    @property
    def descriptor_key(self):
        return self.make_key()

    def calculate(self, mol):
        r = defaultdict(lambda: defaultdict(int))
        for a in mol.GetAtoms():
            if a.GetAtomicNum() != 6:
                continue

            double = 0
            triple = 0
            carbon = 0

            for b in a.GetBonds():
                other = b.GetBeginAtom()
                if a.GetIdx() == other.GetIdx():
                    other = b.GetEndAtom()

                if other.GetAtomicNum() == 6:
                    carbon += 1

                bt = b.GetBondTypeAsDouble()
                if bt == 2.0:
                    double += 1
                elif bt == 3.0:
                    triple += 1

            if (double == 2 and triple == 0) or (triple == 1 and double == 0):
                SP = 1
            elif double == 1 and triple == 0:
                SP = 2
            else:
                SP = 3

            r[SP][carbon] += 1

        return r


class CarbonTypes(CarbonTypesBase):
    descriptor_defaults = [
        (1, 1), (2, 1),
        (1, 2), (2, 2), (3, 2),
        (1, 3), (2, 3), (3, 3), (4, 3),
    ]

    @property
    def descriptor_name(self):
        return 'C{}SP{}'.format(self.nCarbon, self.SP)

    @property
    def descriptor_key(self):
        return self.make_key(self.nCarbon, self.SP)

    @property
    def dependencies(self):
        return dict(CT=CarbonTypesCache.make_key())

    def __init__(self, nCarbon, SP):
        self.nCarbon = nCarbon
        self.SP = SP

    def calculate(self, mol, CT):
        return CT[self.SP][self.nCarbon]


class HybridizationRatio(CarbonTypesBase):
    descriptor_defaults = [()]

    descriptor_name = 'HybRatio'

    @property
    def descriptor_key(self):
        return self.make_key()

    @property
    def dependencies(self):
        return dict(CT=CarbonTypesCache.make_key())

    def calculate(self, mol, CT):
        Nsp3 = float(sum(CT[3].values()))
        Nsp2 = float(sum(CT[2].values()))
        return Nsp3 / (Nsp2 + Nsp3)
