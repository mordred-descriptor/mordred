import numpy as np


class MordredException(Exception):
    critical = False

    def __float__(self):
        return np.nan


class MordredValueError(MordredException, ValueError):
    pass
