import numpy as np

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
    assert isinstance(df.copy(), MordredDataFrame)


def test_fill_missing():
    f = df.fill_missing()
    assert f.iloc[0, 0], 1
    assert f.iloc[0, 1], 2
    assert f.iloc[1, 0], 1
    assert np.isnan(f.iloc[1, 1])
    assert np.isnan(f.iloc[2, 0])
    assert f.iloc[2, 1], 4
