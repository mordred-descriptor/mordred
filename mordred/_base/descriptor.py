from abc import ABCMeta, abstractmethod
from inspect import isabstract
from contextlib import contextmanager

import six
import numpy as np

import operator


class MissingValueException(Exception):
    "internally used exception"

    __slots__ = ('error',)

    def __init__(self, error):
        self.error = error


class Descriptor(six.with_metaclass(ABCMeta, object)):
    r"""abstract base class of descriptors."""

    __slots__ = '_context',

    explicit_hydrogens = True
    kekulize = False
    require_connected = False
    require_3D = False

    def __reduce_ex__(self, version):
        return self.__class__, self.parameters()

    @abstractmethod
    def parameters(self):
        '''get __init__ arguments of this descriptor instance.

        (abstract method)
        '''
        raise NotImplementedError('not implemented Descriptor.parameters method')

    @property
    def as_argument(self):
        '''argument representation of descriptor'''
        return self

    @staticmethod
    def _pretty(v):
        v = getattr(v, 'as_argument', v)
        return repr(v)

    def __repr__(self):
        return '{}({})'.format(
            self.__class__,
            ', '.join(self._pretty(a) for a in self.parameters())
        )

    def __hash__(self):
        return hash((self.__class__, self.parameters()))

    def __compare_by_reduce(meth):
        def compare(self, other):
            l = self.__class__, self.parameters()
            r = other.__class__, other.parameters()
            return getattr(l, meth)(r)

        return compare

    __eq__ = __compare_by_reduce('__eq__')
    __ne__ = __compare_by_reduce('__ne__')

    __lt__ = __compare_by_reduce('__lt__')
    __gt__ = __compare_by_reduce('__gt__')
    __le__ = __compare_by_reduce('__le__')
    __ge__ = __compare_by_reduce('__ge__')

    rtype = None

    @classmethod
    def preset(cls):
        r"""generate preset descriptor instances.

        :rtype: iterable
        """
        return ()

    def dependencies(self):
        r"""descriptor dependencies.

        :rtype: {:py:class:`str`: (:py:class:`Descriptor` or :py:class:`None`)} or :py:class:`None`
        """
        pass

    @abstractmethod
    def calculate(self):
        r"""calculate descriptor value.

        (abstract method)
        """
        raise NotImplementedError('not implemented Descriptor.calculate method')

    @classmethod
    def is_descriptor_class(cls, desc):
        r"""check calculatable descriptor class or not.

        :rtype: :py:class:`bool`
        """
        return (
            isinstance(desc, type) and
            issubclass(desc, cls) and
            not isabstract(desc)
        )

    @property
    def mol(self):
        '''get molecule'''
        return self._context.get_mol(self)

    @property
    def coord(self):
        '''get 3D coordinate'''
        if not self.require_3D:
            self.fail(AttributeError('use 3D coordinate in 2D descriptor'))

        return self._context.get_coord(self)

    @contextmanager
    def rethrow_zerodiv(self):
        '''treat zero div as known exception'''
        with np.errstate(divide='raise', invalid='raise'):
            try:
                yield
            except (FloatingPointError, ZeroDivisionError) as e:
                self.fail(ZeroDivisionError(*e.args))

    def fail(self, exception):
        '''raise known exception and return missing value'''
        raise MissingValueException(exception)

    @contextmanager
    def rethrow_na(self, exception):
        '''treat any exceptions as known exception'''
        try:
            yield
        except exception as e:
            self.fail(e)

    def _unary_common(name, operator):
        def unary(self):
            return UnaryOperatingDescriptor(name.format(self), operator, self)

        return unary

    def _binary_common(name, operator):
        def binary(self, other):
            if not isinstance(other, Descriptor):
                other = ConstDescriptor(str(other), other)

            return BinaryOperatingDescriptor(name.format(self, other), operator, self, other)

        return binary

    __add__ = _binary_common('({}+{})', operator.add)
    __sub__ = _binary_common('({}-{})', operator.sub)
    __mul__ = _binary_common('({}*{})', operator.mul)
    __truediv__ = _binary_common('({}/{})', operator.truediv)
    __floordiv__ = _binary_common('({}//{})', operator.floordiv)
    __mod__ = _binary_common('({}%{})', operator.mod)
    __pow__ = _binary_common('({}**{})', operator.pow)

    __neg__ = _unary_common('-{}', operator.neg)
    __pos__ = _unary_common('+{}', operator.pos)
    __abs__ = _unary_common('|{}|', operator.abs)

    __ceil__ = _unary_common('⌈{}⌉', np.ceil)
    __floor__ = _unary_common('⌊{}⌋', np.floor)
    __trunc__ = _unary_common('trunc({})', np.trunc)


class UnaryOperatingDescriptor(Descriptor):
    @classmethod
    def preset(cls):
        return cls()

    def parameters(self):
        return self.name, self.operator, self.value

    def __init__(self, name, operator, value):
        self.name = name
        self.operator = operator
        self.value = value

    def __str__(self):
        return self.name

    def dependencies(self):
        return {
            'value': self.value,
        }

    def calculate(self, value):
        return self.operator(value)


class ConstDescriptor(Descriptor):
    @classmethod
    def preset(cls):
        return cls()

    def parameters(self):
        return self.name, self.value

    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __str__(self):
        return self.name

    def calculate(self):
        return self.value


class BinaryOperatingDescriptor(Descriptor):
    @classmethod
    def preset(cls):
        return cls()

    def parameters(self):
        return self.name, self.operator, self.left, self.right

    def __init__(self, name, operator, left, right):
        self.name = name
        self.operator = operator
        self.left = left
        self.right = right

    def __str__(self):
        return self.name

    def dependencies(self):
        return {
            'left': self.left,
            'right': self.right,
        }

    def calculate(self, left, right):
        return self.operator(left, right)
