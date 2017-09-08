import numpy as np
from pandas import Series, DataFrame

from .util import is_missing

__all__ = ("MordredDataFrame", "Series")


class MordredDataFrame(DataFrame):
    @property
    def _constructor(self):
        return self.__class__

    def fill_missing(self, value=np.nan, inplace=False):
        t = self if inplace else self.copy()

        t[t.applymap(is_missing)] = value
        return t
