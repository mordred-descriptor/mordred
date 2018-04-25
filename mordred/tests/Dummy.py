from mordred import Descriptor


class Dummy1(Descriptor):
    pass


class Dummy2(Dummy1):
    since = "1.0.0"

    def parameters(self):
        return ()

    def calculate(self):
        return 0.0

    rtype = float


class Dummy3(Dummy2):
    @classmethod
    def preset(cls, version):
        yield cls()


class Dummy4(Dummy3):
    def __init__(self, i=0):
        self.i = i

    def __eq__(self, other):
        return isinstance(other, Dummy4) and self.i == other.i

    def __str__(self):
        return "Dummy4_{}".format(self.i)

    @classmethod
    def preset(cls, version):
        yield cls(0)
        yield cls(1)
