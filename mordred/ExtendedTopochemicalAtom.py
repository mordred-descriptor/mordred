r"""Extended Topochemical Atom(ETA) descriptor.

References
    * :doi:`10.1021/ci0342066`
    * :doi:`10.1016/j.jhazmat.2013.03.023`

"""
from __future__ import division

import numpy as np
from rdkit import Chem

from . import _atomic_property as ap
from ._base import Descriptor
from ._util import atoms_to_numpy
from .RingCount import RingCount

__all__ = (
    "EtaCoreCount",
    "EtaShapeIndex",
    "EtaVEMCount",
    "EtaCompositeIndex",
    "EtaFunctionalityIndex",
    "EtaBranchingIndex",
    "EtaDeltaAlpha",
    "EtaEpsilon",
    "EtaDeltaEpsilon",
    "EtaDeltaBeta",
    "EtaPsi",
    "EtaDeltaPsi",
)


class AlterMolecule(Descriptor):
    __slots__ = ("explicit_hydrogens", "_saturated")

    kekulize = True
    require_connected = True

    def parameters(self):
        return self.explicit_hydrogens, self._saturated

    def __init__(self, explicit_hydrogens, saturated=False):
        self._saturated = saturated
        self.explicit_hydrogens = explicit_hydrogens

    def __str__(self):
        b = "Saturated" if self._saturated else "Reference"
        H = "H" if self.explicit_hydrogens else ""

        return "{}Mol{}".format(b, H)

    def calculate(self):
        new = Chem.RWMol(Chem.Mol())
        ids = {}
        for a in self.mol.GetAtoms():
            if a.GetAtomicNum() == 1:
                continue

            if self._saturated:
                new_a = Chem.Atom(a.GetAtomicNum())
                new_a.SetFormalCharge(a.GetFormalCharge())

            else:
                new_a = Chem.Atom(6)

            ids[a.GetIdx()] = new.AddAtom(new_a)

        for bond in self.mol.GetBonds():
            ai = bond.GetBeginAtom()
            aj = bond.GetEndAtom()

            if not self._saturated and (ai.GetDegree() > 4 or aj.GetDegree() > 4):
                self.fail(ValueError("bond degree greater then 4"))

            i = ids.get(ai.GetIdx())
            j = ids.get(aj.GetIdx())

            if i is not None and j is not None:
                if self._saturated and (
                    ai.GetAtomicNum() != 6 or aj.GetAtomicNum() != 6
                ):
                    order = bond.GetBondType()
                else:
                    order = Chem.BondType.SINGLE

                new.AddBond(i, j, order)

        new = Chem.Mol(new)
        if Chem.SanitizeMol(new, catchErrors=True) != 0:
            typ = "saturated" if self._saturated else "referense"
            self.fail(ValueError("cannot sanitize {} mol".format(typ)))

        if self.explicit_hydrogens:
            new = Chem.AddHs(new)

        Chem.Kekulize(new)

        return new


class EtaBase(Descriptor):
    __slots__ = ()

    explicit_hydrogens = False
    kekulize = True
    require_connected = True

    rtype = float


class EtaCoreCount(EtaBase):
    r"""ETA core count descriptor.

    .. math::

        \alpha = \sum_{i = 1}^A \frac{Z_i - Z_i^v}{Z_i^v} \cdot \frac{1}{PN_i - 1}

    where :math:`Z_i` and :math:`Z_i^v` are number of total and valence electons,
    :math:`PN` is periodic number.

    :type averaged: bool
    :param averaged: averaged by number of heavy count

    :type reference: bool
    :param reference: use reference alkane
        (same graph structure, but all atoms are carbon and all bonds are single bond)

    :returns: reference and valence of any atoms > 4
    """

    since = "1.0.0"
    __slots__ = ("_averaged", "_reference")

    def description(self):
        return "{}ETA core count{}".format(
            "averaged " if self._averaged else "",
            " for reference graph" if self._reference else "",
        )

    @classmethod
    def preset(cls, version):
        return map(cls, [False, True])

    def __str__(self):
        suffix = "_R" if self._reference else ""
        ave = "A" if self._averaged else ""

        return "{}ETA_alpha{}".format(ave, suffix)

    def parameters(self):
        return self._averaged, self._reference

    def __init__(self, averaged=False, reference=False):
        self._averaged = averaged
        self._reference = reference

    def dependencies(self):
        if self._reference:
            return {"rmol": AlterMolecule(self.explicit_hydrogens)}

    def calculate(self, rmol=None):
        mol = rmol if self._reference else self.mol

        v = sum(ap.get_core_count(a) for a in mol.GetAtoms())
        if self._averaged:
            v /= mol.GetNumAtoms()

        return v


class EtaShapeIndex(EtaBase):
    r"""ETA shape index descriptor.

    .. math::
        {\rm shape}_t = \frac{\alpha_t}{\alpha}

    where :math:`\alpha_t` is p(alpha value only atoms which bond to 1 heavy atom),
    y(3), or x(4).

    :type type: str
    :param type: one of shape_types
    """

    since = "1.0.0"
    __slots__ = ("_type",)

    def description(self):
        return "ETA shape index (type: {})".format(self._type)

    shape_types = ("p", "y", "x")
    _type_to_degree = {"p": 1, "y": 3, "x": 4}

    @classmethod
    def preset(cls, version):
        return (cls(t) for t in cls.shape_types)

    def __str__(self):
        return "ETA_shape_{}".format(self._type)

    def parameters(self):
        return (self._type,)

    def __init__(self, type="p"):
        assert type in self.shape_types

        self._type = type

    def dependencies(self):
        return {"a": EtaCoreCount(False)}

    def calculate(self, a):
        d = self._type_to_degree[self._type]

        return (
            sum(ap.get_core_count(a) for a in self.mol.GetAtoms() if a.GetDegree() == d)
            / a
        )


class EtaVEMCount(EtaBase):
    r"""ETA VEM(valence electron mobile) count descriptor.

    .. math::
        \beta^{\rm s} = \frac{1}{2} \sum^A_{i=1} \beta^{\rm s}_i

        \beta^{\rm s}_i = \sum^A_{j = 1} x_{ij}\sigma_{ij}

        x_{ij} = \begin{cases}
            0.5  & \left( \left| \epsilon_i - \epsilon_j \right| \leq 0.3 \right) \\
            0.75 & \left( \left| \epsilon_i - \epsilon_j \right| > 0.3 \right)
        \end{cases}

        \epsilon_i = - \alpha_i + 0.3 Z^{\rm v}

    where :math:`\sigma_{ij}` is sigma bond count between i and j.

    .. math::
        \beta^{\rm ns\delta} = \sum^A_{i = 1} \beta^{\rm ns\delta}_i

    where :math:`\beta^{\rm ns\delta}_i` is
    0.5 if i-th atom is making resonance with an aromatic ring.

    .. math::
        \beta^{\rm ns} = \frac{1}{2} \sum^A_{i=1} \beta^{\rm ns}_i

        \beta^{\rm ns}_i = \sum^A_{j = 1} y_{ij}\pi_{ij} + \beta^{\rm ns\delta}_i

        y_{ij} = \begin{cases}
            2.0 & \left( {\rm ij\ is\ aromatic\ bond} \right) \\
            1.5 & \left( \left| \epsilon_i - \epsilon_j \right| > 0.3 \right) \\
            1.0 & \left( \left| \epsilon_i - \epsilon_j \right| \leq 0.3 \right)
        \end{cases}

    where :math:`\pi_{ij}` is pi bond count between i and j.

    .. math::
        \beta = \beta^{\rm s} + \beta^{\rm ns}

    :type type: str
    :param type: one of beta_types

    :type averaged: bool
    :param averaged: averaged by heavy atom count
    """

    since = "1.0.0"
    __slots__ = ("_type", "_averaged")

    def description(self):
        if self._type == "":
            d = ""
        elif self._type == "s":
            d = "sigma contribution to "
        elif self._type == "ns":
            d = "nonsigma contribution to "
        else:
            d = "delta contribution to "

        return "{}{}valence electron mobile count".format(
            "averaged " if self._averaged else "", d
        )

    @classmethod
    def preset(cls, version):
        return (cls(b, a) for b in cls.beta_types for a in [False, True])

    def __str__(self):
        typ = "_{}".format(self._type) if self._type else ""
        ave = "A" if self._averaged else ""

        return "{}ETA_beta{}".format(ave, typ)

    beta_types = ("", "s", "ns", "ns_d")

    def parameters(self):
        return self._type, self._averaged

    def __init__(self, type="", averaged=False):
        assert type in self.beta_types
        self._type = type
        self._averaged = averaged

    def _get_beta_s(self, atom):
        return ap.get_eta_beta_sigma(atom) / 2.0

    def _get_beta_ns_d(self, atom):
        v = ap.get_eta_beta_delta(atom)
        return v

    def _get_beta_ns(self, atom):
        return ap.get_eta_beta_non_sigma(atom) / 2.0 + self._get_beta_ns_d(atom)

    def _get_beta_(self, atom):
        return self._get_beta_s(atom) + self._get_beta_ns(atom)

    def calculate(self):
        getter = getattr(self, "_get_beta_" + self._type)

        if getter:
            v = sum(getter(a) for a in self.mol.GetAtoms())

        if self._averaged:
            v /= self.mol.GetNumAtoms()

        return v


class EtaCompositeIndex(EtaBase):
    r"""ETA composite index descriptor.

    .. math::
        \eta = \sum_{i < j} \left( \frac{\gamma_i \gamma_j}{r_{ij}^2} \right)^{0.5}

        \gamma_i = \frac{\alpha_i}{\beta_i}

    where :math:`r_{ij}` is graph distance.

    .. math::
        \eta^{\rm local} = \sum_{i < j, r_{ij} = 1} \left( \gamma_i \gamma_j \right)^{0.5}

    :type reference: bool
    :param reference: use reference alkane.

    :type local: bool
    :param local: use :math:`\eta^{\rm local}`

    :type averaged: bool
    :param averaged: averaged

    :returns: reference and valence of any atoms > 4
    """

    since = "1.0.0"
    __slots__ = ("_reference", "_local", "_averaged")

    def description(self):
        return "{}{}ETA composite index{}".format(
            "averaged " if self._averaged else "",
            "local " if self._local else "",
            " for reference graph",
        )

    @classmethod
    def preset(cls, version):
        ft = [False, True]
        return (cls(r, l, a) for r in ft for l in ft for a in ft)

    def __str__(self):
        suffixes = []

        if self._reference:
            suffixes.append("R")

        if self._local:
            suffixes.append("L")

        if len(suffixes) > 0:
            suffix = "_" + "".join(suffixes)
        else:
            suffix = ""

        ave = "A" if self._averaged else ""

        return "{}ETA_eta{}".format(ave, suffix)

    def parameters(self):
        return self._reference, self._local, self._averaged

    def __init__(self, reference=False, local=False, averaged=False):
        self._reference = reference
        self._local = local
        self._averaged = averaged

    def dependencies(self):
        if self._reference:
            return {"rmol": AlterMolecule(self.explicit_hydrogens)}

    def calculate(self, rmol=None):
        mol = rmol if self._reference else self.mol
        D = Chem.GetDistanceMatrix(mol, force=True)

        if self._local:

            def checker(r):
                return r == 1

        else:

            def checker(r):
                return r != 0

        gamma = atoms_to_numpy(ap.get_eta_gamma, mol)

        v = sum(
            sum(
                np.sqrt(gamma[i] * gamma[j] / r ** 2)
                for j, r in enumerate(row)
                if i < j and checker(r)
            )
            for i, row in enumerate(D)
        )

        if self._averaged:
            v /= mol.GetNumAtoms()

        return v


class EtaFunctionalityIndex(EtaBase):
    r"""ETA functionality index descriptor.

    .. math::
        \eta^{\rm F} = \eta^{\rm R} - \eta

    where :math:`\eta^{\rm R}` is eta value for reference alkane.

    :type local: bool
    :param local: use local eta

    :type averaged: bool
    :param averaged: averaged
    """

    since = "1.0.0"
    __slots__ = ("_local", "_averaged")

    def description(self):
        return "{}{}ETA functionality index".format(
            "averaged " if self._averaged else "", "local " if self._local else ""
        )

    @classmethod
    def preset(cls, version):
        return (cls(l, a) for l in [False, True] for a in [False, True])

    def __str__(self):
        loc = "L" if self._local else ""
        ave = "A" if self._averaged else ""

        return "{}ETA_eta_F{}".format(ave, loc)

    def parameters(self):
        return self._local, self._averaged

    def __init__(self, local=False, averaged=False):
        self._local = local
        self._averaged = averaged

    def dependencies(self):
        return {
            "eta": EtaCompositeIndex(local=self._local),
            "eta_R": EtaCompositeIndex(local=self._local, reference=True),
        }

    def calculate(self, eta, eta_R):
        v = eta_R - eta
        if self._averaged:
            v /= self.mol.GetNumAtoms()

        return v


class EtaBranchingIndex(EtaBase):
    r"""ETA branching index descriptor.

    .. math::
        \eta^{\rm B} = \eta^{\rm local,N} - \eta^{local,R} + 0.086 N^{\rm R}

    where :math:`\eta^{\rm local,N}` is :math:`\eta^{\rm local}` for normal alkane.
    :math:`N^{\rm R}` is ring count.

    :type ring: bool
    :param ring: use ring count or not

    :type averaged: bool
    :param averaged: averaged

    :returns: NaN when A < 2
    """

    since = "1.0.0"
    __slots__ = ("_ring", "_averaged")

    def description(self):
        return "{}ETA branching index{}".format(
            "averaged " if self._averaged else "",
            " (use ring count)" if self._ring else "",
        )

    @classmethod
    def preset(cls, version):
        return (cls(r, a) for r in [False, True] for a in [False, True])

    def __str__(self):
        ring = "R" if self._ring else ""
        ave = "A" if self._averaged else ""

        return "{}ETA_eta_B{}".format(ave, ring)

    def parameters(self):
        return self._ring, self._averaged

    def __init__(self, ring=True, averaged=False):
        self._ring = ring
        self._averaged = averaged

    def dependencies(self):
        return {
            "NR": RingCount() if self._ring else None,
            "eta_RL": EtaCompositeIndex(reference=True, local=True),
        }

    def calculate(self, eta_RL, NR):
        N = self.mol.GetNumAtoms()

        if N <= 1:
            self.fail(ValueError("single atom"))
        elif N == 2:
            eta_NL = 1.0
        else:
            eta_NL = np.sqrt(2) + 0.5 * (N - 3)

        v = eta_NL - eta_RL + 0.086 * (NR or 0)

        if self._averaged:
            v /= N

        return v


class EtaDeltaAlpha(EtaBase):
    r"""ETA delta alpha descriptor.

    .. math::
        \Delta\alpha_{\rm A} = \max\left(\frac{\alpha - \alpha^{\rm R}}{A}, 0\right)

        \Delta\alpha_{\rm B} = \max\left(\frac{\alpha^{\rm R} - \alpha}{A}, 0\right)

    :type type: str
    :param type: one of delta_types
    """

    since = "1.0.0"
    __slots__ = ("_type",)

    def description(self):
        return "ETA delta alpha (type: {})".format(self._type)

    delta_types = ("A", "B")

    @classmethod
    def preset(cls, version):
        return (cls(t) for t in cls.delta_types)

    def __str__(self):
        return "ETA_dAlpha_{}".format(self._type)

    def parameters(self):
        return (self._type,)

    def __init__(self, type="A"):
        assert type in self.delta_types

        self._type = type

    def dependencies(self):
        return {"alpha": EtaCoreCount(), "alpha_R": EtaCoreCount(reference=True)}

    def calculate(self, alpha, alpha_R):
        if self._type == "A":
            d = alpha - alpha_R
        else:
            d = alpha_R - alpha

        return max(d / self.mol.GetNumAtoms(), 0.0)


class EtaEpsilon(EtaBase):
    r"""ETA epsilon descriptor.

    .. math::
        \epsilon^i = \frac{\epsilon^i}{N^i} (i \leq 4)
        \epsilon^5 = \frac{\epsilon^2 + \epsilon^{\rm XH}}{N^2 + N^{\rm XH}}

    types(i)
        1
            all atoms
        2
            heavy atoms
        3
            all atoms of reference alkane
        4
            all atoms of saturated carbon skeleton(reduce C-C bonds)
        XH
            hydrogens bond to hetero atoms

    :type type: str
    :param type: one of epsilon_types

    :returns: type = 3 and valence of any atoms > 4
    """

    since = "1.0.0"
    __slots__ = ("_type",)

    def description(self):
        return "ETA epsilon (type: {})".format(self._type)

    @classmethod
    def preset(cls, version):
        return map(cls, cls.epsilon_types)

    def __str__(self):
        return "ETA_epsilon_{}".format(self._type)

    @property
    def explicit_hydrogens(self):
        return self._type != 2

    epsilon_types = tuple(range(1, 6))

    def parameters(self):
        return (self._type,)

    def __init__(self, type=1):
        self._type = type

    def dependencies(self):
        if self._type == 3:
            return {"rmol": AlterMolecule(self.explicit_hydrogens)}
        elif self._type == 4:
            return {"rmol": AlterMolecule(self.explicit_hydrogens, True)}

    def calculate(self, rmol=None):
        mol = rmol if self._type in [3, 4] else self.mol

        if self._type == 5:
            eps = [
                ap.get_eta_epsilon(a)
                for a in mol.GetAtoms()
                if a.GetAtomicNum() != 1 or a.GetNeighbors()[0].GetAtomicNum() != 6
            ]
            return sum(eps) / len(eps)

        return sum(ap.get_eta_epsilon(a) for a in mol.GetAtoms()) / mol.GetNumAtoms()


class EtaDeltaEpsilon(EtaBase):
    r"""ETA delta epsilon descriptor.

    .. math::
        \Delta \epsilon^{\rm A} = \epsilon^1 - \epsilon^3

        \Delta \epsilon^{\rm B} = \epsilon^1 - \epsilon^4

        \Delta \epsilon^{\rm C} = \epsilon^3 - \epsilon^4

        \Delta \epsilon^{\rm D} = \epsilon^2 - \epsilon^5

    :type type: str
    :param type: one of delta_epsilon_types
    """

    since = "1.0.0"
    __slots__ = ("_type",)

    def description(self):
        return "ETA delta epsilon (type: {})".format(self._type)

    @classmethod
    def preset(cls, version):
        return map(cls, cls.delta_epsilon_types)

    def __str__(self):
        return "ETA_dEpsilon_{}".format(self._type)

    delta_epsilon_types = tuple("ABCD")

    def parameters(self):
        return (self._type,)

    def __init__(self, type="A"):
        self._type = type

    _deps = {"A": (1, 3), "B": (1, 4), "C": (3, 4), "D": (2, 5)}

    def dependencies(self):
        L, R = self._deps[self._type]
        return {"L": EtaEpsilon(L), "R": EtaEpsilon(R)}

    def calculate(self, L, R):
        return L - R


class EtaDeltaBeta(EtaBase):
    r"""ETA delta beta descriptor.

    .. math::
        \Delta\beta = \beta^{\rm ns} - \beta^{\rm s}

    :type averaged: bool
    :param averaged: averaged
    """

    since = "1.0.0"
    __slots__ = ("_averaged",)

    def description(self):
        return "{}ETA delta beta".format("averaged " if self._averaged else "")

    @classmethod
    def preset(cls, version):
        return (cls(a) for a in [False, True])

    def __str__(self):
        ave = "A" if self._averaged else ""

        return "{}ETA_dBeta".format(ave)

    def parameters(self):
        return (self._averaged,)

    def __init__(self, averaged=False):
        self._averaged = averaged

    def dependencies(self):
        return {"ns": EtaVEMCount("ns"), "s": EtaVEMCount("s")}

    def calculate(self, ns, s):
        v = ns - s

        if self._averaged:
            v /= self.mol.GetNumAtoms()

        return v


class EtaPsi(EtaBase):
    r"""ETA psi descriptor.

    .. math::
        \psi_1 = \frac{\alpha}{A \cdot \epsilon^2}
    """

    since = "1.0.0"
    __slots__ = ()

    def description(self):
        return "ETA psi"

    @classmethod
    def preset(cls, version):
        yield cls()

    def __str__(self):
        return "ETA_psi_1"

    def parameters(self):
        return ()

    def dependencies(self):
        return {"a": EtaCoreCount(), "e": EtaEpsilon(2)}

    def calculate(self, a, e):
        return a / (self.mol.GetNumAtoms() * e)


class EtaDeltaPsi(EtaBase):
    r"""ETA delta psi descriptor.

    .. math::
        \Delta\psi_{\rm A} = \max\left(0.714 - \psi_1, 0\right)

        \Delta\psi_{\rm B} = \max\left(\psi_1 - 0.714, 0\right)

    :type type: str
    :param type: one of delta_psi_types
    """

    since = "1.0.0"
    __slots__ = ("_type",)

    def description(self):
        return "ETA delta psi (type: {})".format(self._type)

    @classmethod
    def preset(cls, version):
        return map(cls, cls.delta_psi_types)

    def __str__(self):
        return "ETA_dPsi_{}".format(self._type)

    delta_psi_types = ("A", "B")

    def parameters(self):
        return (self._type,)

    def __init__(self, type="A"):
        assert type in self.delta_psi_types

        self._type = type

    def dependencies(self):
        return {"psi": EtaPsi()}

    def calculate(self, psi):
        L = 0.714
        R = psi

        if self._type == "B":
            L, R = R, L

        return max(L - R, 0.0)
