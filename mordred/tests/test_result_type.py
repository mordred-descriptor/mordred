try:
    from collections.abc import Iterable
except ImportError:
    import warnings

    warnings.filterwarnings(
        "ignore",
        ".* ABCs from 'collections' .*",
        DeprecationWarning,
    )
    from collections import Iterable

import numpy as np

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


result = Result(None, [1, error.Error(ValueError("error"), [])], [Dummy1(), Dummy2()])


def test_length():
    assert len(result) == 2


def test_fill_missing():
    assert np.isnan(result.fill_missing()[1])


def test_drop_missing():
    assert len(result.drop_missing()) == 1


def test_items():
    i = result.items()
    assert isinstance(i, Iterable)
    li = list(i)
    assert len(li) == 2
    assert len(li[0]) == 2

    d, v = li[0]
    assert isinstance(d, Descriptor)
    assert isinstance(v, int)


def test_keys():
    i = result.keys()
    assert isinstance(i, Iterable)
    li = list(i)
    assert len(li), 2

    assert isinstance(li[0], Descriptor)


def test_iter():
    i = iter(result)
    assert isinstance(i, Iterable)
    li = list(i)
    assert len(li) == 2

    assert isinstance(li[0], int)


def test_reversed():
    i = reversed(result)
    assert isinstance(i, Iterable)
    li = list(i)
    assert len(li) == 2

    assert isinstance(li[1], int)


def test_asdict_non_rawkey():
    d = result.asdict()
    assert len(d) == 2
    assert isinstance(d, dict)
    for k in d.keys():
        assert isinstance(k, str)


def test_asdict_rawkey():
    d = result.asdict(True)
    assert len(d) == 2
    assert isinstance(d, dict)
    for k in d.keys():
        assert isinstance(k, Descriptor)


def test_ix():
    assert 1 == result.ix[0]
    assert is_missing(result.ix[1])


def test_name():
    assert 1 == result.name["Dummy1"]
    assert is_missing(result.name["Dummy2"])


def test_getitem():
    assert 1 == result[0]
    assert is_missing(result[1])
    assert 1 == result["Dummy1"]
    assert is_missing(result["Dummy2"])
