from collections import Iterator

from mordred import (
    Calculator,
    descriptors,
    get_descriptors_in_module,
    get_descriptors_from_module,
)
from nose.tools import ok_


def test_descriptor_order():
    calc = Calculator(descriptors)
    it = iter(calc.descriptors)
    before = next(it).__module__
    for current in it:
        current = current.__module__
        assert before <= current, "{!r} > {!r}".format(before, current)

        before = current


def test_get_descriptors_in_module():
    old = get_descriptors_from_module(descriptors, True)
    new = get_descriptors_in_module(descriptors, True)

    yield isinstance, new, Iterator
    yield ok_, len(list(new)) > 100

    for a, b in zip(old, new):
        assert a == b
