from mordred import Calculator, Descriptor
from nose.tools import eq_, ok_, raises

from . import Dummy


def test_Calculator_descriptors():
    calc = Calculator()

    def check(l, msg=None):
        yield eq_, len(calc), l, msg
        yield ok_, all(isinstance(d, Descriptor) for d in calc.descriptors)

    # register instance
    calc.register(Dummy.Dummy2())
    for c in check(1, "instance register failed"):
        yield c

    # register class
    yield raises(ValueError)(lambda: calc.register(Dummy.Dummy1))
    for c in check(1):
        yield c

    calc.register(Dummy.Dummy2)
    for c in check(1):
        yield c

    calc.register(Dummy.Dummy3)
    for c in check(2):
        yield c

    calc.register(Dummy.Dummy4)
    for c in check(4):
        yield c

    # delete
    del calc.descriptors
    for c in check(0):
        yield c

    # register module
    calc.register(Dummy)
    for c in check(3):
        yield c

    eq_(calc["Dummy4_1"], Dummy.Dummy4(1))

    # set instance
    calc.descriptors = Dummy.Dummy2()
    for c in check(1):
        yield c
