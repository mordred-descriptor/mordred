from mordred import Descriptor


class Dummy1(Descriptor):
    pass


class Dummy2(Dummy1):
    def as_key(self):
        return self.__class__, ()

    def calculate(self, mol):
        return 0.0

    rtype = float


class Dummy3(Dummy2):
    @classmethod
    def preset(cls):
        yield cls()


class Dummy4(Dummy3):
    @classmethod
    def preset(cls):
        yield cls()
        yield cls()
