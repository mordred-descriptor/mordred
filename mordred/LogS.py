# flake8: noqa: S1
# -*- coding: UTF-8 -*-

from __future__ import unicode_literals

import math

from rdkit import Chem

from ._base import Descriptor
from .Weight import Weight

__all__ = ("LogS",)


_smarts_logs = {
    "[NH0;X3;v3]": 0.71535,
    "[NH2;X3;v3]": 0.41056,
    "[nH0;X3]": 0.82535,
    "[OH0;X2;v2]": 0.31464,
    "[OH0;X1;v2]": 0.14787,
    "[OH1;X2;v2]": 0.62998,
    "[CH2;!R]": -0.35634,
    "[CH3;!R]": -0.33888,
    "[CH0;R]": -0.21912,
    "[CH2;R]": -0.23057,
    "[ch0]": -0.37570,
    "[ch1]": -0.22435,
    "F": -0.21728,
    "Cl": -0.49721,
    "Br": -0.57982,
    "I": -0.51547,
}


_smarts_logs_molecules = [
    (Chem.MolFromSmarts(smarts), log) for smarts, log in _smarts_logs.items()
]


class LogS(Descriptor):
    r"""Filter-it™ LogS descriptor.

    http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/filter-it/1.0.2/filter-it.html#installation
    """

    __slots__ = ()

    since = "1.1.0"
    explicit_hydrogens = False
    kekulize = False

    @classmethod
    def preset(cls, version):
        yield cls()

    def dependencies(self):
        return {"MW": Weight(exact=False)}

    def description(self):
        return "Filter-it™ LogS"

    def __str__(self):
        return "FilterItLogS"

    def parameters(self):
        return ()

    def calculate(self, MW):
        logS = 0.89823 - 0.10369 * math.sqrt(MW)
        for smarts, log in _smarts_logs_molecules:
            logS += len(self.mol.GetSubstructMatches(smarts)) * log

        return logS

    rtype = float
