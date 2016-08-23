import numpy as np
import six
from abc import ABCMeta, abstractproperty


class MissingValueBase(six.with_metaclass(ABCMeta, object)):
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

    @abstractproperty
    def header(self):
        raise NotImplementedError('require header')


class Missing(MissingValueBase):
    __slots__ = ()
    header = 'Missing'


class Error(MissingValueBase):
    __slots__ = ()
    header = 'ERROR'


class MordredException(Exception):
    __slots__ = ()


class MultipleFragments(MordredException):
    __slots__ = ()

    def __str__(self):
        return 'multiple fragments'


class Missing3DCoordinate(MordredException):
    __slots__ = ()

    def __str__(self):
        return 'missing 3D coordinate'
