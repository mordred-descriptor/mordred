from rdkit import Chem

from .._util import conformer_to_numpy
from ..error import Missing3DCoordinate


class Context(object):
    __slots__ = "_mols", "_coords", "n_frags", "name", "_stack"

    def __init__(self, mols, coords, n_frags, name):
        self._mols = mols
        self._coords = coords
        self.n_frags = n_frags
        self.name = name

    def __reduce_ex__(self, version):
        return self.__class__, (self._mols, self._coords, self.n_frags, self.name)

    def __str__(self):
        return self.name

    __tf = {True, False}

    @classmethod
    def from_query(cls, mol, require_3D, explicit_hydrogens, kekulizes, id):
        if not isinstance(mol, Chem.Mol):
            raise TypeError("{!r} is not rdkit.Chem.Mol instance".format(mol))

        n_frags = len(Chem.GetMolFrags(mol))

        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        else:
            name = Chem.MolToSmiles(Chem.RemoveHs(mol))

        mols, coords = {}, {}

        for eh, ke in ((eh, ke) for eh in explicit_hydrogens for ke in kekulizes):
            m = Chem.AddHs(mol) if eh else Chem.RemoveHs(mol)

            if ke:
                Chem.Kekulize(m)

            if require_3D:
                try:
                    conf = m.GetConformer(id)
                    if conf.Is3D():
                        coords[eh, ke] = conformer_to_numpy(conf)
                except ValueError:
                    pass

            m.RemoveAllConformers()
            mols[eh, ke] = m

        return cls(mols, coords, n_frags, name)

    @classmethod
    def from_calculator(cls, calc, mol, id):
        return cls.from_query(mol, calc._require_3D, calc._explicit_hydrogens, calc._kekulizes, id)

    def get_coord(self, desc):
        try:
            return self._coords[desc.explicit_hydrogens, desc.kekulize]
        except KeyError:
            desc.fail(Missing3DCoordinate())

    def get_mol(self, desc):
        return self._mols[desc.explicit_hydrogens, desc.kekulize]

    def reset(self):
        self._stack = []

    def add_stack(self, d):
        self._stack.append(d)

    def get_stack(self):
        return self._stack
