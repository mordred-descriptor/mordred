from ._base import Descriptor


class FragmentComplexity(Descriptor):
    r'''
    fragment complexity descriptor

    .. math::
        {\rm fragCpx} = \left| B^2 - A^2 + A \right| + \frac{H}{100}

    where
    :math:`A` is number of atoms,
    :math:`B` is number of bonds,
    :math:`H` is number of hetero atoms

    Returns:
        float: fragment complexity
    '''

    explicit_hydrogens = False
    descriptor_name = 'fragCpx'

    def calculate(self, mol):
        A = mol.GetNumAtoms()
        B = mol.GetNumBonds()
        H = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 6)
        return abs(B ** 2 - A ** 2 + A) + float(H) / 100
