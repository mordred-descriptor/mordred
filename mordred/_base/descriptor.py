import six
from abc import ABCMeta, abstractmethod
from inspect import isabstract


class Descriptor(six.with_metaclass(ABCMeta, object)):
    r"""abstract base class of descriptors."""

    explicit_hydrogens = True
    kekulize = False
    require_connected = False
    require_3D = False

    _reduce_ex_version = 3

    @abstractmethod
    def __reduce_ex__(self, version):
        pass

    @property
    def as_argument(self):
        return self

    @staticmethod
    def _pretty(v):
        v = getattr(v, 'as_argument', v)
        return repr(v)

    def __repr__(self):
        cls, args = self.__reduce_ex__(self._reduce_ex_version)
        return '{}({})'.format(cls.__name__, ', '.join(self._pretty(a) for a in args))

    def __hash__(self):
        return hash(self.__reduce_ex__(self._reduce_ex_version))

    def __compare_by_reduce(meth):
        def compare(self, other):
            l = self.__reduce_ex__(self._reduce_ex_version)
            r = other.__reduce_ex__(self._reduce_ex_version)
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

    # def __call__(self, mol, coord_id=-1):
    #     r"""calculate single descriptor value.

    #     :returns: descriptor result
    #     :rtype: scalar
    #     """
    #     return Calculator(self)(mol, coord_id)[0]

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
