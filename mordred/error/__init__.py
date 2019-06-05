"""Error objects."""

from abc import ABCMeta, abstractproperty

import six
import numpy as np


class MissingValueBase(six.with_metaclass(ABCMeta, object)):
    """Base class of missing values.

    Args:
        error (Exception): error object
        stack (callstack)

    """

    __slots__ = "error", "stack"

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
        return "{} ({})".format(self.error, "/".join(str(d) for d in self.stack))

    @abstractproperty
    def header(self):
        """Header of warning message.

        Returns:
            str

        """
        raise NotImplementedError("require header")


class Missing(MissingValueBase):
    """known errored value."""

    __slots__ = ()
    header = "Missing"


class Error(MissingValueBase):
    """unknown errored value."""

    __slots__ = ()
    header = "ERROR"


class MordredException(Exception):
    __slots__ = ()


class MultipleFragments(MordredException):
    """multiple fragments detected on require_connected Descriptor."""

    __slots__ = ()

    def __str__(self):
        return "multiple fragments"


class Missing3DCoordinate(MordredException):
    """missing 3D coordinate on require_3D Descriptor."""

    __slots__ = ()

    def __str__(self):
        return "missing 3D coordinate"


class DuplicatedDescriptorName(MordredException):
    """duplicated string replisantation of descriptor."""

    __slots__ = ()

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __str__(self):
        return "duplicated descriptor name: {!r} and {!r}".format(self.a, self.b)


class Timeout(MordredException):
    """calculation timed out."""

    __slots__ = ()

    def __str__(self):
        return "timed out"
