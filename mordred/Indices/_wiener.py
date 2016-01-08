from .._base import Descriptor
from .._common import DistanceMatrix


class Wiener(Descriptor):
    explicit_hydrogens = False

    descriptor_defaults = [(False,), (True,)]

    @property
    def dependencies(self):
        return dict(
            D=DistanceMatrix.make_key(
                self.explicit_hydrogens,
                False,
                False,
            )
        )

    @property
    def descriptor_name(self):
        return 'WPol' if self.polarity else 'WPath'

    def __init__(self, polarity=False):
        self.polarity = polarity

    @property
    def descriptor_key(self):
        return self.make_key(self.polarity)

    def calculate(self, mol, D):
        if self.polarity:
            return int(0.5 * (D == 3).sum())
        else:
            return int(0.5 * D.sum())
