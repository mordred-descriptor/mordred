import numpy as np

from nose.tools import eq_, ok_
from mordred.error import Error, Missing
from mordred._base.pandas_module import MordredDataFrame

df = MordredDataFrame(
    [
        [1, 2],
        [1, Error(ValueError("error"), ())],
        [Missing(ValueError("missing"), ()), 4],
    ]
)


def test_copy():
    ok_(isinstance(df.copy(), MordredDataFrame))


def test_fill_missing():
    f = df.fill_missing()
    yield eq_, f.iloc[0, 0], 1
    yield eq_, f.iloc[0, 1], 2
    yield eq_, f.iloc[1, 0], 1
    yield ok_, np.isnan(f.iloc[1, 1])
    yield ok_, np.isnan(f.iloc[2, 0])
    yield eq_, f.iloc[2, 1], 4
