from rdkit import Chem
from .._util import conformer_to_numpy


class Context(object):
    __slots__ = '_mols', '_coords', 'n_frags', 'name', '_stack'

    def __setstate__(self, dict):
        self._mols = dict.get('_mols', {})
        self._coords = dict.get('_coords', {})
        self.n_frags = dict.get('n_frags', -1)
        self.name = dict.get('name', 'unknown')

    def __reduce_ex__(self, version):
        return self.__class__, (None,), {
            '_mols': self._mols,
            '_coords': self._coords,
            'n_frags': self.n_frags,
            'name': self.name,
        }

    def __str__(self):
        return self.name

    @classmethod
    def from_calculator(cls, calc, mol, id):
        return cls(mol, calc._require_3D, calc._explicit_hydrogens, calc._kekulizes, id)

    __tf = set([True, False])

    def __init__(self, mol, require_3D=False, explicit_hydrogens=__tf, kekulizes=__tf, id=-1):
        if mol is None:
            return

        self._mols = {}
        self._coords = {}

        self.n_frags = len(Chem.GetMolFrags(mol))

        if mol.HasProp('_Name'):
            self.name = mol.GetProp('_Name')
        else:
            self.name = Chem.MolToSmiles(Chem.RemoveHs(mol))

        for eh, ke in ((eh, ke) for eh in explicit_hydrogens for ke in kekulizes):
            m = (Chem.AddHs if eh else Chem.RemoveHs)(mol)

            if ke:
                Chem.Kekulize(m)

            if require_3D:
                self._coords[eh, ke] = conformer_to_numpy(m.GetConformer(id))

            m.RemoveAllConformers()
            self._mols[eh, ke] = m

    def get_coord(self, explicit_hydrogens, kekulize):
        return self._coords[explicit_hydrogens, kekulize]

    def get_mol(self, explicit_hydrogens, kekulize):
        return self._mols[explicit_hydrogens, kekulize]

    def reset(self):
        self._stack = []

    def add_stack(self, d):
        self._stack.append(d)

    def get_stack(self):
        return self._stack
