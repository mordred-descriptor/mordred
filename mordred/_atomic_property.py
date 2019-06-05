from __future__ import division

import os

import six
import numpy as np
from rdkit import Chem
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges

from ._base import Descriptor
from ._util import atoms_to_numpy


halogen = {9, 17, 35, 53, 85, 117}


getter_list = []

if six.PY2:
    from collections import OrderedDict

    getters = OrderedDict()
else:
    getters = {}


def getter(short, **attrs):
    def proc(f):
        f.short = short

        for a, v in attrs.items():
            setattr(f, a, v)

        getters[short] = f
        getter_list.append(f)
        return f

    if short in getters:
        raise ValueError("duplicated short name of atomic property")

    return proc


@getter(short="c", long="gasteiger charge", gasteiger_charges=True)
def get_gasteiger_charge(atom):
    return (
        atom.GetDoubleProp("_GasteigerCharge") + atom.GetDoubleProp("_GasteigerHCharge")
        if atom.HasProp("_GasteigerHCharge")
        else 0.0
    )


class PeriodicTable(object):
    __slots__ = ("data",)

    _datadir = os.path.join(os.path.dirname(__file__), "data")

    def __init__(self, data=None):
        self.data = data

    @classmethod
    def load(cls, name, conv=float):
        def read(v):
            if "-" in v:
                return np.nan

            try:
                return conv(v)
            except ValueError:
                return

        self = cls()

        with open(os.path.join(cls._datadir, name)) as file:
            self.data = [
                v for v in (read(l.split("#")[0]) for l in file) if v is not None
            ]

        return self

    def __getitem__(self, i):
        if i < 1:
            return np.nan

        try:
            return self.data[i - 1]
        except IndexError:
            return np.nan

    def map(self, f):
        new = self.__class__()
        new.data = [f(d) for d in self.data]
        return new


mass = PeriodicTable.load("mass.txt")
vdw_radii = PeriodicTable.load("van_der_waals_radii.txt")
vdw_volume = vdw_radii.map(lambda r: 4.0 / 3.0 * np.pi * r ** 3)
sanderson = PeriodicTable.load("sanderson_electron_negativity.txt")
pauling = PeriodicTable.load("pauling_electron_negativity.txt")
allred_rocow = PeriodicTable.load("allred_rocow_electron_negativity.txt")
polarizability94 = PeriodicTable.load("polarizalibity94.txt")
polarizability78 = PeriodicTable.load("polarizalibity78.txt")
ionization_potentials = PeriodicTable.load("ionization_potential.txt")
period = PeriodicTable(
    ([1] * 2)
    + ([2] * 8)
    + ([3] * 8)
    + ([4] * 18)
    + ([5] * 18)
    + ([6] * 32)
    + ([7] * 32)
)

mc_gowan_volume = PeriodicTable.load("mc_gowan_volume.txt")


_table = Chem.GetPeriodicTable()


def GetElementSymbol(i):
    return _table.GetElementSymbol(i)


if six.PY2:

    def GetAtomicNumber(symbol):
        if isinstance(symbol, unicode):  # noqa: F821
            symbol = str(symbol)
        return _table.GetAtomicNumber(symbol)


else:

    def GetAtomicNumber(symbol):
        return _table.GetAtomicNumber(symbol)


# http://dx.doi.org/10.1002%2Fjps.2600721016
@getter(short="dv", long="valence electrons", valence=True)
def get_valence_electrons(atom):
    N = atom.GetAtomicNum()
    if N == 1:
        return 0

    Zv = _table.GetNOuterElecs(N) - atom.GetFormalCharge()
    Z = atom.GetAtomicNum() - atom.GetFormalCharge()
    hi = atom.GetTotalNumHs()
    he = sum(1 for a in atom.GetNeighbors() if a.GetAtomicNum() == 1)
    h = hi + he
    return (Zv - h) / (Z - Zv - 1)


@getter(short="d", long="sigma electrons", valence=True)
def get_sigma_electrons(atom):
    return sum(1 for a in atom.GetNeighbors() if a.GetAtomicNum() != 1)


# http://www.edusoft-lc.com/molconn/manuals/400/chaptwo.html
# p. 283
@getter(short="s", long="intrinsic state", require_connected=True, valence=True)
def get_intrinsic_state(atom):
    i = atom.GetAtomicNum()
    d = get_sigma_electrons(atom)
    dv = get_valence_electrons(atom)

    if d == 0:
        return np.nan

    return ((2.0 / period[i]) ** 2 * dv + 1) / d


def get_core_count(atom):
    Z = atom.GetAtomicNum()
    if Z == 1:
        return 0.0

    Zv = _table.GetNOuterElecs(Z)
    PN = period[Z]

    return (Z - Zv) / (Zv * (PN - 1))


def get_eta_epsilon(atom):
    Zv = _table.GetNOuterElecs(atom.GetAtomicNum())
    return 0.3 * Zv - get_core_count(atom)


def get_eta_beta_sigma(atom):
    e = get_eta_epsilon(atom)
    return sum(
        0.5 if abs(get_eta_epsilon(a) - e) <= 0.3 else 0.75
        for a in atom.GetNeighbors()
        if a.GetAtomicNum() != 1
    )


def get_eta_nonsigma_contribute(bond):
    if bond.GetBondType() is Chem.BondType.SINGLE:
        return 0.0

    f = 1.0
    if bond.GetBondTypeAsDouble() == Chem.BondType.TRIPLE:
        f = 2.0

    a = bond.GetBeginAtom()
    b = bond.GetEndAtom()

    dEps = abs(get_eta_epsilon(a) - get_eta_epsilon(b))

    if bond.GetIsAromatic():
        y = 2.0
    elif dEps > 0.3:
        y = 1.5
    else:
        y = 1.0

    return y * f


def get_eta_beta_delta(atom):
    if (
        atom.GetIsAromatic()
        or atom.IsInRing()
        or _table.GetNOuterElecs(atom.GetAtomicNum()) - atom.GetTotalValence() <= 0
    ):
        return 0.0

    for b in atom.GetNeighbors():
        if b.GetIsAromatic():
            return 0.5

    return 0.0


def get_other_atom(bond, atom):
    begin = bond.GetBeginAtom()
    if atom.GetIdx() != begin.GetIdx():
        return begin

    return bond.GetEndAtom()


def get_eta_beta_non_sigma(atom):
    return sum(
        get_eta_nonsigma_contribute(b)
        for b in atom.GetBonds()
        if get_other_atom(b, atom).GetAtomicNum() != 1
    )


def get_eta_gamma(atom):
    beta = (
        get_eta_beta_sigma(atom)
        + get_eta_beta_non_sigma(atom)
        + get_eta_beta_delta(atom)
    )
    if beta == 0:
        return np.nan

    return get_core_count(atom) / beta


@getter(short="Z", long="atomic number")
def get_atomic_number(a):
    return a.GetAtomicNum()


@getter(short="m", long="mass")
def get_mass(a):
    return mass[a.GetAtomicNum()]


@getter(short="v", long="vdw volume")
def get_vdw_volume(a):
    return vdw_volume[a.GetAtomicNum()]


@getter(short="se", long="sanderson EN")
def get_sanderson_en(a):
    return sanderson[a.GetAtomicNum()]


@getter(short="pe", long="pauling EN")
def get_pauling_en(a):
    return pauling[a.GetAtomicNum()]


@getter(short="are", long="allred-rocow EN")
def get_allred_rocow_en(a):
    return allred_rocow[a.GetAtomicNum()]


@getter(short="p", long="polarizability")
def get_polarizability(a):
    return polarizability94[a.GetAtomicNum()]


@getter(short="i", long="ionization potential")
def get_ionization_potential(a):
    return ionization_potentials[a.GetAtomicNum()]


def get_mc_gowan_volume(a):
    return mc_gowan_volume[a.GetAtomicNum()]


def get_properties(charge=False, valence=False):
    for f in getters.values():
        if not charge and getattr(f, "gasteiger_charges", False):
            continue

        if not valence and getattr(f, "valence", False):
            continue

        yield f.short


class AtomicProperty(Descriptor):
    __slots__ = "explicit_hydrogens", "prop", "_initialized"

    def __str__(self):
        return "Prop{}".format(self.as_argument)

    def get_long(self):
        return getattr(self.prop, "long", self.prop.__name__)

    @property
    def as_argument(self):
        return getattr(self.prop, "short", self.prop.__name__)

    def parameters(self):
        return self.explicit_hydrogens, self.prop

    def __new__(cls, explicit_hydrogens, prop):
        if isinstance(prop, cls):
            prop._initialized = True
            return prop

        return super(AtomicProperty, cls).__new__(cls)

    def __init__(self, explicit_hydrogens, prop):
        if getattr(self, "_initialized", False):
            return

        self.explicit_hydrogens = explicit_hydrogens
        self.prop = getters.get(prop)

        if self.prop is not None:
            return

        if callable(prop):
            self.prop = prop
            return

        raise TypeError("atomic property is not callable: {!r}".format(prop))

    def calculate(self):
        if getattr(self.prop, "gasteiger_charges", False):
            ComputeGasteigerCharges(self.mol)

        r = atoms_to_numpy(self.prop, self.mol)

        nans = np.isnan(r)
        if np.any(nans):
            atms = set(np.array([a.GetSymbol() for a in self.mol.GetAtoms()])[nans])
            self.fail(
                ValueError("missing {} for {}".format(self.get_long(), list(atms)))
            )

        return r

    _carbon = Chem.Atom(6)

    @property
    def carbon(self):
        return self.prop(self._carbon)
