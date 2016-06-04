import six
from abc import ABCMeta, abstractmethod
from inspect import isabstract


class Descriptor(six.with_metaclass(ABCMeta, object)):
    r"""abstract base class of descriptors."""

    __slots__ = '_context',

    explicit_hydrogens = True
    kekulize = False
    require_connected = False
    require_3D = False

    def __reduce_ex__(self, version):
        return self.as_key()

    @abstractmethod
    def as_key(self):
        raise TypeError('not implemented Descriptor.as_key method')

    @property
    def as_argument(self):
        return self

    @staticmethod
    def _pretty(v):
        v = getattr(v, 'as_argument', v)
        return repr(v)

    def __repr__(self):
        cls, args = self.as_key()
        return '{}({})'.format(cls.__name__, ', '.join(self._pretty(a) for a in args))

    def __hash__(self):
        return hash(self.as_key())

    def __compare_by_reduce(meth):
        def compare(self, other):
            l = self.as_key()
            r = other.as_key()
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
    def calculate(self, mol):
        r"""calculate descriptor value.

        (abstract method)
        """
        raise TypeError('not implemented Descriptor.calculate method')

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
        return self._context.get_mol(self.explicit_hydrogens, self.kekulize)

    @property
    def coord(self):
        return self._context.get_coord(self.explicit_hydrogens, self.kekulize)
