import numpy as np


class Error(object):
    __slots__ = 'error', 'stack'

    def __reduce_ex__(self, version):
        return self.__class__, (self.error, self.stack)

    def __init__(self, error, stack):
        self.error = error
        self.stack = stack

    def __float__(self):
        return np.nan

    def __add__(self, other):
        return np.nan

    def __sub__(self, other):
        return np.nan

    def __str__(self):
        return '{} ({})'.format(self.error, '/'.join(str(d) for d in self.stack))


class MordredException(Exception):
    __slots__ = ()
    pass


class MultipleFragments(MordredException):
    __slots__ = ()

    def __init__(self):
        pass

    def __str__(self):
        return 'multiple fragments'
