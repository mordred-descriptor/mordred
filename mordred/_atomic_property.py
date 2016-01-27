from numpy import nan, pi

from rdkit import Chem


def attr(**attrs):
    def proc(f):
        for a, v in attrs.items():
            setattr(f, a, v)

        return f

    return proc


@attr(gasteiger_charges=True, name='c')
def get_charge_explicitHs(atom):
    return atom.GetDoubleProp('_GasteigerCharge')


@attr(gasteiger_charges=True, name='c')
def get_charge_implicitHs(atom):
    return atom.GetDoubleProp('_GasteigerCharge') +\
        atom.GetDoubleProp('_GasteigerHCharge')


na = nan

# Handbook of Chemistry and Physics, 94th Edition, 2013-2014, pg1-11
mass = [
    nan,
    # 1
    1.008, 4.002602,
    # 2
    6.94, 9.012182, 10.81, 12.011, 14.007, 15.999, 18.9984032, 20.1797,
    # 3
    22.98976928, 24.3050, 26.9815386, 28.085, 30.973762, 32.06, 35.45, 39.948,
    # 4
    39.0983, 40.078, 44.955912, 47.867, 50.9415, 51.9961, 54.938045, 55.845, 58.933195,
    58.6934, 63.546, 65.38, 69.723, 72.63, 74.92160, 78.96, 79.904, 83.798,
    # 5
    85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.96, 98, 101.07, 102.90550,
    106.42, 107.8682, 112.411, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.293,
    # 6
    132.9054519, 137.327,
    138.90547, 140.116, 140.90765, 144.242, 145, 150.36, 151.964, 157.25,
    158.92535, 162.500, 164.93032, 167.259, 168.93421, 173.054, 174.9668,
    178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084,
    196.966569, 200.59, 204.38, 207.2, 208.98040, 210, 210, 222,
    # 7
    223, 226,
    227, 232.03806, 231.03588, 238.02891, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 262,
    261, 262, 266, 264, 269, 268, 271
]


# Handbook of Chemistry and Physics, 94th Edition, 2013-2014, pg9-49
Rvdw = [
    nan,
    # 1
    1.10, 1.40,
    # 2
    1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
    # 3
    2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88,
    # 4
    2.75, 2.31, 2.15, 2.11, 2.07, 2.06, 2.05, 2.04, 2.00,
    1.97, 1.96, 2.01, 1.87, 2.11, 1.85, 1.90, 1.85, 2.02,
    # 5
    3.03, 2.49, 2.32, 2.23, 2.18, 2.17, 2.16, 2.13, 2.10,
    2.10, 2.11, 2.18, 1.93, 2.17, 2.06, 2.06, 1.98, 2.16,
    # 6
    3.43, 2.68,
    2.43, 2.42, 2.40, 2.39, 2.38, 2.36, 2.35, 2.34, 2.33, 2.31, 2.30, 2.29, 2.27, 2.26, 2.24,
    2.23, 2.22, 2.18, 2.16, 2.16, 2.13, 2.13, 2.14, 2.23, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20,
    # 7
    3.48, 2.83,
    2.47, 2.45, 2.43, 2.41, 2.39, 2.43, 2.44, 2.45, 2.44, 2.45, 2.45, 2.45, 2.46, 2.46, 2.46,
]

Vvdw = [4. / 3. * pi * r ** 3 for r in Rvdw]

# Sanderson, RT. Electronegativity and Bond Energy.  J. Am. Chem. Soc. 1983, 105: 2259-2261.
# Remaining values come from http://www.talete.mi.it/help/dragon_help/weighting_schemes.htm
Sanderson = [
    nan,
    # 1
    2.592, na,
    # 2
    0.670, 1.810, 2.275, 2.746, 3.194, 3.654, 4.000, 4.5,
    # 3
    0.560, 1.318, 1.714, 2.138, 2.515, 2.957, 3.475, 3.31,
    # 4
    0.445, 0.946, 1.02, 1.09, 1.39, 1.66, 2.2, 2.2, 2.56,
    1.94, 2.033, 2.223, 2.419, 2.618, 2.816, 3.014, 3.219, 2.91,
    # 5
    0.312, 0.721, 0.65, 0.9, 1.42, 1.15, na, na, na,
    na, 1.826, 1.978, 2.138, 2.298, 2.458, 2.618, 2.778, 2.34,
    # 6
    0.220, 0.651,
    na, na, na, na, na, na, na, na, na, na, na, na, na, na, na,
    na, na, 0.98, na, na, na, na, na, 2.195, 2.246, 2.291, 2.342
]

# https://github.com/cdk/cdk/blob/master/misc/extra/src/main/resources/org/openscience/cdk/config/data/electroneg-pauling.txt

Pauling = [
    nan,
    # 1
    2.2, na,
    # 2
    0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, na,
    # 3
    0.93, 1.31, 1.61, 1.9, 2.19, 2.58, 3.16, na,
    # 4
    0.82, 1.0, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88,
    1.91, 1.9, 1.65, 1.81, 2.01, 2.18, 2.55, 2.96, 3.0,
    # 5
    0.82, 0.95, 1.22, 1.33, 1.6, 2.16, 1.9, 2.2, 2.28,
    2.2, 1.93, 1.69, 1.78, 1.96, 2.05, 2.1, 2.66, 2.6,
    # 6
    0.79, 0.89,
    1.1, 1.12, 1.13, 1.14, na, 1.17, na, 1.2, na, 1.22, 1.23, 1.24, 1.25, na, 1.27,
    1.3, 1.5, 2.36, 1.9, 2.2, 2.2, 2.28, 2.54, 2.0, 1.62, 2.33, 2.02, 2.0, 2.2, na,

    # 7
    0.7, 0.9,
    1.1, 1.3, 1.5, 1.38, 1.36, 1.28, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3,
]

# http://www.chm.davidson.edu/ronutt/che115/electroneg.htm
Allred_Rocow = [
    nan,
    # 1
    2.20, na,
    # 2
    0.97, 1.47, 2.01, 2.50, 3.07, 3.50, 4.10, na,
    # 3
    1.01, 1.23, 1.47, 1.74, 2.06, 2.44, 2.83, na,
    # 4
    0.91, 1.04, 1.20, 1.32, 1.45, 1.56, 1.60, 1.64, 1.70,
    1.75, 1.75, 1.66, 1.82, 2.02, 2.20, 2.48, 2.74, na,
    # 5
    0.89, 0.99, 1.11, 1.22, 1.23, 1.30, 1.36, 1.42, 1.45,
    1.35, 1.42, 1.46, 1.49, 1.72, 1.82, 2.01, 2.21, na,
    # 6
    0.86, 0.97,
    1.08, na, na, na, na, na, na, na, na, na, na, na, na, na, na,
    1.23, 1.33, 1.40, 1.46, 1.52, 1.55, 1.44, 1.42, 1.44, 1.44, 1.55, 1.67, 1.76, 1.90
]


# Handbook of Chemistry and Physics, 94th Edition, 2013-2014, pg10-188
Polarizabilities94 = [
    nan,
    # 1
    0.666793, 0.2050522,
    # 2
    24.33, 5.60, 3.03, 1.67, 1.10, 0.802, 0.557, 0.39432,
    # 3
    24.11, 10.6, 6.8, 5.53, 3.63, 2.90, 2.18, 1.6411,
    # 4
    43.06, 22.8, 17.8, 14.6, 12.4, 11.6, 9.4, 8.4, 7.5,
    6.8, 6.2, 5.75, 8.12, 5.84, 4.31, 3.77, 3.05, 2.4844,
    # 5
    47.24, 23.5, 22.7, 17.9, 15.7, 12.8, 11.4, 9.6, 8.6,
    4.8, 6.78, 7.36, 10.2, 7.84, 6.6, 5.5, 5.35, 4.044,
    # 6
    59.42, 39.7,
    31.1, 29.6, 28.2, 31.4, 30.1, 28.8, 27.7, 23.5, 25.5, 24.5, 23.6, 22.7, 21.8, 20.9, 21.9,
    16.2, 13.1, 11.1, 9.7, 8.5, 7.6, 6.5, 5.8, 5.02, 7.6, 7.01, 7.4, 6.8, 6.0, 5.3,
    # 7
    48.6, 38.3,
    32.1, 32.1, 25.4, 24.9, 24.8, 24.5, 23.3, 23.0, 22.7, 20.5, 19.7, 23.8, 18.2, 16.4, na,
    na, na, na, na, na, na, na, na, 4.06, na, 4.59, na, na, na, na,
    # 8
    24.26,
]

# Handbook of Chemistry and Physics, 78th Edition
Polarizabilities78 = [
    nan,
    # 1
    0.666793, 0.204956,
    # 2
    24.3, 5.6, 3.03, 1.76, 1.1, 0.802, 0.557, 0.3956,
    # 3
    23.6, 10.6, 6.8, 5.38, 3.63, 2.9, 2.18, 1.6411,

    # 4
    43.4, 22.8, 17.8, 14.6, 12.4, 11.6, 9.4, 8.4, 7.5,
    6.8, 6.1, 7.1, 8.12, 6.07, 4.31, 3.77, 3.05, 2.4844,

    # 5
    47.3, 27.6, 22.7, 17.9, 15.7, 12.8, 11.4, 9.6, 8.6,
    4.8, 7.2, 7.2, 10.2, 7.7, 6.6, 5.5, 5.35, 4.044,

    # 6
    59.6, 39.7,
    31.1, 29.6, 28.2, 31.4, 30.1, 28.8, 27.7, 23.5, 25.5, 24.5, 23.6, 22.7, 21.8, 21.0, 21.9,
    16.2, 13.1, 11.1, 9.7, 8.5, 7.6, 6.5, 5.8, 5.7, 7.6, 6.8, 7.4, 6.8, 6.0, 5.3,

    # 7
    48.7, 38.3,
    32.1, 32.1, 25.4, 27.4, 24.8, 24.5, 23.3, 23.0, 22.7, 20.5, 19.7, 23.8, 18.2, 17.5
]


# Handbook of Chemistry and Physics, 94th Edition, 2013-2014, pg10-197
Ionpotentials = [
    nan,
    # 1
    13.598443, 24.587387,
    # 2
    5.391719, 9.32270, 8.29802, 11.26030, 14.5341, 13.61805, 17.4228, 21.56454,
    # 3
    5.139076, 7.646235, 5.985768, 8.15168, 10.48669, 10.36001, 12.96763, 15.759610,
    # 4
    4.3406633, 6.11316, 6.56149, 6.82812, 6.74619, 6.76651, 7.43402, 7.9024, 7.88101,
    7.6398, 7.72638, 9.394199, 5.999301, 7.89943, 9.7886, 9.75239, 11.8138, 13.99961,
    # 5
    4.177128, 5.69485, 6.2173, 6.63390, 6.75885, 7.09243, 7.28, 7.36050, 7.45890,
    8.3369, 7.57623, 8.99382, 5.78636, 7.34392, 8.60839, 9.0096, 10.45126, 12.12984,
    # 6
    3.893905, 5.211664,
    5.5769, 5.5387, 5.473, 5.5250, 5.582, 5.6437, 5.67038, 6.14980,
    5.8638, 5.9389, 6.0215, 6.1077, 6.18431, 6.25416, 5.42586,
    6.82507, 7.54957, 7.86403, 7.83352, 8.43823, 8.96702, 8.9588,
    9.22553, 10.4375, 6.108194, 7.41663, 7.2855, 8.414, na, 10.7485,
    # 7
    4.072741, 5.278423,
    5.17, 6.3067, 5.89, 6.1941, 6.2657, 6.0260, 5.9738,
    5.9914, 6.1979, 6.2817, 6.42, 6.50, 6.58, 6.65, 4.9,
    6.0
]

period = [
    nan,
    1, 1,
    2, 2, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 3, 3, 3, 3, 3,

    4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4,

    5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5,

    6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,

    7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
]


McGowanVolume = [
    nan,
    # 1
    8.71, 6.75,
    # 2
    22.23, 20.27, 18.31, 16.35, 14.39, 12.43, 10.47, 8.51,
    # 3
    32.71, 30.75, 28.79, 26.83, 24.87, 22.91, 20.95, 18.99,
    # 4
    51.89, 50.28, 48.68, 47.07, 45.47, 43.86, 42.26, 40.65, 39.05,
    37.44, 35.84, 34.23, 32.63, 31.02, 29.42, 27.81, 26.21, 24.60,
    # 5
    60.22, 58.61, 57.01, 55.40, 53.80, 52.19, 50.59, 48.98, 47.38,
    45.77, 44.17, 42.56, 40.96, 39.35, 37.75, 36.14, 34.54, 32.93,
    # 6
    77.25, 76.00,
    74.75, 73.49, 72.24, 70.99, 69.74, 68.49, 67.23, 65.98,
    64.73, 63.48, 62.23, 60.97, 59.72, 58.47, 57.22,
    55.97, 54.71, 53.46, 52.21, 50.96, 49.71, 48.45, 47.20,
    45.95, 44.70, 43.45, 42.19, 40.94, 39.69, 38.44,
    # 7
    75.59, 74.34,
    73.09, 71.83, 70.58, 69.33, 68.08, 66.83, 65.57, 64.32,
    63.07, 61.82, 60.57, 59.31, 58.06, 56.81, 55.56,
]


table = Chem.GetPeriodicTable()


# http://dx.doi.org/10.1002%2Fjps.2600721016
@attr(name='delta_v')
def get_valence_electrons(atom):
    N = atom.GetAtomicNum()
    if N == 1:
        return 0

    Zv = table.GetNOuterElecs(N) - atom.GetFormalCharge()
    Z = atom.GetAtomicNum() - atom.GetFormalCharge()
    hi = atom.GetTotalNumHs()
    he = sum(1 for a in atom.GetNeighbors() if a.GetAtomicNum() == 1)
    h = hi + he
    return float(Zv - h) / float(Z - Zv - 1)


@attr(name='delta')
def get_sigma_electrons(atom):
    return sum(1 for a in atom.GetNeighbors()
               if a.GetAtomicNum() != 1)


# http://www.edusoft-lc.com/molconn/manuals/400/chaptwo.html
# p. 283
@attr(require_connected=True, name='s')
def get_intrinsic_state(atom):
    i = atom.GetAtomicNum()
    d = get_sigma_electrons(atom)
    dv = get_valence_electrons(atom)

    if d == 0:
        return nan

    return ((2. / period[i]) ** 2 * dv + 1) / d


def get_core_count(atom):
    Z = atom.GetAtomicNum()
    if Z == 1:
        return 0.0

    Zv = table.GetNOuterElecs(Z)
    PN = period[Z]

    return float(Z - Zv) / (Zv * (PN - 1))


def get_eta_epsilon(atom):
    Zv = table.GetNOuterElecs(atom.GetAtomicNum())
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
    if atom.GetIsAromatic() or\
            atom.IsInRing() or\
            table.GetNOuterElecs(atom.GetAtomicNum()) - atom.GetTotalValence() <= 0:
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
    beta = get_eta_beta_sigma(atom) + get_eta_beta_non_sigma(atom) + get_eta_beta_delta(atom)
    if beta == 0:
        return nan

    return get_core_count(atom) / beta


@attr(name='Z')
def get_atomic_number(a):
    return a.GetAtomicNum()


@attr(name='m')
def get_mass(a):
    return mass[a.GetAtomicNum()]


@attr(name='v')
def get_vdw_volume(a):
    return Vvdw[a.GetAtomicNum()]


@attr(name='e')
def get_sanderson_en(a):
    return Sanderson[a.GetAtomicNum()]


@attr(name='pe')
def get_pauling_en(a):
    return Pauling[a.GetAtomicNum()]


@attr(name='are')
def get_allred_rocow_en(a):
    return Allred_Rocow[a.GetAtomicNum()]


@attr(name='p')
def get_polarizability(a):
    return Polarizabilities94[a.GetAtomicNum()]


@attr(name='i')
def get_ionpotential(a):
    return Ionpotentials[a.GetAtomicNum()]


def get_mc_gowan_volume(a):
    return McGowanVolume[a.GetAtomicNum()]


getters = dict(
    Z=get_atomic_number,
    m=get_mass,
    v=get_vdw_volume,
    e=get_sanderson_en,
    se=get_sanderson_en,
    pe=get_pauling_en,
    are=get_allred_rocow_en,
    p=get_polarizability,
    i=get_ionpotential,
    s=get_intrinsic_state,
    delta=get_sigma_electrons,
    delta_v=get_valence_electrons,
)


def get_properties(charge=False, istate=False):
    if charge:
        yield 'c'

    for p in ['Z', 'm', 'v', 'e', 'pe', 'are', 'p', 'i']:
        yield p

    if istate:
        yield 's'


def getter(p, explicit_hydrogens):
    if p in getters:
        p = getters[p]

    if p == 'c':
        if explicit_hydrogens:
            p = get_charge_explicitHs
        else:
            p = get_charge_implicitHs

    if hasattr(p, '__call__'):
        return getattr(p, 'name', p.__name__), p

    raise ValueError('atomic property is not callable: {!r}'.format(p))
