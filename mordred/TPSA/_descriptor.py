from .._base import Descriptor
from rdkit.Chem import rdMolDescriptors, BondType
from collections import defaultdict


class TPSA(Descriptor):
    @classmethod
    def preset(cls):
        return [(True,), (False,)]

    @property
    def descriptor_name(self):
        return 'TPSA(NO)' if self.no_only else 'TPSA'

    def __init__(self, no_only=True):
        self.no_only = no_only

    @property
    def descriptor_key(self):
        return self.make_key(self.no_only)

    def calculate(self, mol):
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        if self.no_only:
            return tpsa

        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()

            if atomic_num == 15:
                tpsa += self._get_phosphorus_contrib(atom)
            elif atomic_num == 16:
                tpsa += self._get_sulfur_contrib(atom)

        return tpsa

    @staticmethod
    def _hydrogen_count(atom):
        return atom.GetTotalNumHs() + \
            sum(1 for a in atom.GetNeighbors() if a.GetAtomicNum() == 1)

    @staticmethod
    def _bond_type_count(atom):
        cnt = defaultdict(int)
        for bond in atom.GetBonds():
            if bond.GetIsAromatic():
                cnt[BondType.AROMATIC] += 1
            else:
                cnt[bond.GetBondType()] += 1

        return dict(cnt)

    @classmethod
    def _get_phosphorus_contrib(self, atom):
        nH = self._hydrogen_count(atom)
        cnt = self._bond_type_count(atom)

        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            return 0.0

        if nH == 1 and cnt == {BondType.SINGLE: 3, BondType.DOUBLE: 1}:
            return 23.47
        elif nH == 0:
            if cnt == {BondType.SINGLE: 3}:
                return 13.59
            elif cnt == {BondType.SINGLE: 1, BondType.DOUBLE: 1}:
                return 34.14
            elif cnt == {BondType.SINGLE: 3, BondType.DOUBLE: 1}:
                return 9.81

        return 0.0

    @classmethod
    def _get_sulfur_contrib(self, atom):
        nH = self._hydrogen_count(atom)
        cnt = self._bond_type_count(atom)

        if atom.GetFormalCharge() != 0:
            return 0.0

        if atom.GetIsAromatic():
            if nH == 0:
                if cnt == {BondType.AROMATIC: 2}:
                    return 28.24
                elif cnt == {BondType.AROMATIC: 2, BondType.DOUBLE: 1}:
                    return 21.70

        else:
            if nH == 1 and cnt == {BondType.SINGLE: 2}:
                return 38.80
            elif nH == 0:
                if cnt == {BondType.SINGLE: 2}:
                    return 25.30
                elif cnt == {BondType.DOUBLE: 1}:
                    return 32.09
                elif cnt == {BondType.SINGLE: 2, BondType.DOUBLE: 1}:
                    return 19.21
                elif cnt == {BondType.SINGLE: 2, BondType.DOUBLE: 2}:
                    return 8.38

        return 0.0


_descriptors = [TPSA]
__all__ = [d.__name__ for d in _descriptors]
