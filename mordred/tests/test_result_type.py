from collections import Iterable

import numpy as np
from nose.tools import eq_, ok_, raises

from mordred import Result, Descriptor, error, is_missing


class Dummy1(Descriptor):
    def __str__(self):
        return "Dummy1"

    def parameters(self):
        return ()

    def calculate(self):
        return 1


class Dummy2(Descriptor):
    def __str__(self):
        return "Dummy2"

    def parameters(self):
        return ()

    def calculate(self):
        raise ValueError("error")


result = Result([1, error.Error(ValueError("error"), [])], [Dummy1(), Dummy2()])


def test_length():
    eq_(len(result), 2)


def test_fill_missing():
    assert np.isnan(result.fill_missing()[1])


def test_drop_missing():
    eq_(len(result.drop_missing()), 1)


def test_items():
    i = result.items()
    yield ok_, isinstance(i, Iterable)
    li = list(i)
    yield eq_, len(li), 2
    yield eq_, len(li[0]), 2

    d, v = li[0]
    yield ok_, isinstance(d, Descriptor)
    yield ok_, isinstance(v, int)


def test_keys():
    i = result.keys()
    yield ok_, isinstance(i, Iterable)
    li = list(i)
    yield eq_, len(li), 2

    yield ok_, isinstance(li[0], Descriptor)


def test_iter():
    i = iter(result)
    yield ok_, isinstance(i, Iterable)
    li = list(i)
    yield eq_, len(li), 2

    yield ok_, isinstance(li[0], int)


def test_reversed():
    i = reversed(result)
    yield ok_, isinstance(i, Iterable)
    li = list(i)
    yield eq_, len(li), 2

    yield ok_, isinstance(li[1], int)


def test_asdict_non_rawkey():
    d = result.asdict()
    yield eq_, len(d), 2
    yield ok_, isinstance(d, dict)
    for k in d.keys():
        yield ok_, isinstance(k, str)


def test_asdict_rawkey():
    d = result.asdict(True)
    yield eq_, len(d), 2
    yield ok_, isinstance(d, dict)
    for k in d.keys():
        yield ok_, isinstance(k, Descriptor)


def test_ix():
    yield eq_, 1, result.ix[0]
    yield ok_, is_missing(result.ix[1])


def test_name():
    yield eq_, 1, result.name["Dummy1"]
    yield ok_, is_missing(result.name["Dummy2"])


def test_getitem():
    yield eq_, 1, result[0]
    yield ok_, is_missing(result[1])
    yield eq_, 1, result["Dummy1"]
    yield ok_, is_missing(result["Dummy2"])


@raises(TypeError)
def test_getitem_raise():
    result[1.2]
