import inspect
import operator
from abc import ABCMeta, abstractmethod
from contextlib import contextmanager
from packaging.version import Version as StrictVersion

import six
import numpy as np
from rdkit import Chem

if hasattr(inspect, "getfullargspec"):

    def getargs(func):
        return tuple(inspect.getfullargspec(func).args[1:])

else:

    def getargs(func):
        try:
            return tuple(inspect.getargspec(func).args[1:])
        except TypeError:
            return ()


class MissingValueException(Exception):
    """Internally used exception."""

    __slots__ = ("error",)

    def __init__(self, error):
        self.error = error


class DescriptorMeta(ABCMeta):
    def __new__(cls, classname, bases, dict):
        __init__ = dict.get("__init__")
        if __init__ is None:
            for base in bases:
                __init__ = getattr(base, "__init__", None)
                if __init__ is not None:
                    break

        dict["parameter_names"] = getargs(__init__)
        if "since" in dict:
            dict["since"] = StrictVersion(dict["since"])

        return ABCMeta.__new__(cls, classname, bases, dict)


class Descriptor(six.with_metaclass(DescriptorMeta, object)):
    r"""Abstract base class of descriptors.

    Attributes:
        mol(rdkit.Mol): target molecule

    """

    __slots__ = ("_context",)

    explicit_hydrogens = True
    kekulize = False
    require_connected = False
    require_3D = False

    def __reduce_ex__(self, version):
        return self.__class__, self.parameters()

    def description(self):
        pass

    @classmethod
    def preset(cls, version):
        r"""Generate preset descriptor instances.

        Returns:
            Iterable[Descriptor]: preset descriptors

        """
        return ()

    @abstractmethod
    def parameters(self):
        """[abstractmethod] get __init__ arguments of this descriptor instance.

        this method used in pickling and identifying descriptor instance.

        Returns:
            tuple: tuple of __init__ arguments

        """
        raise NotImplementedError("not implemented Descriptor.parameters method")

    def get_parameter_dict(self):
        return dict(zip(self.parameter_names, self.parameters()))

    def to_json(self):
        """Convert to json serializable dictionary.

        Returns:
            dict: dictionary of descriptor

        """
        d, ps = self._to_json()
        if len(ps) == 0:
            return {"name": d}
        else:
            return {"name": d, "args": ps}

    def _to_json(self):
        d = self.__class__.__name__
        ps = self.get_parameter_dict()

        return d, {k: getattr(v, "as_argument", v) for k, v in ps.items()}

    @abstractmethod
    def calculate(self):
        r"""[abstractmethod] calculate descriptor value.

        Returns:
            rtype

        """
        raise NotImplementedError("not implemented Descriptor.calculate method")

    def dependencies(self):
        r"""Descriptor dependencies.

        Returns:
            dict[str, Descriptor or None] or None

        """
        pass

    @property
    def as_argument(self):
        """Argument representation of descriptor.

        Returns:
            any

        """
        return self

    @staticmethod
    def _pretty(v):
        v = getattr(v, "as_argument", v)
        return repr(v)

    def __repr__(self):
        return "{}.{}({})".format(
            self.__class__.__module__,
            self.__class__.__name__,
            ", ".join(self._pretty(a) for a in self.parameters()),
        )

    def __hash__(self):
        return hash((self.__class__, self.parameters()))

    def __compare_by_reduce(meth):
        def compare(self, other):
            L = self.__class__, self.parameters()
            r = other.__class__, other.parameters()
            return getattr(L, meth)(r)

        return compare

    __eq__ = __compare_by_reduce("__eq__")
    __ne__ = __compare_by_reduce("__ne__")

    __lt__ = __compare_by_reduce("__lt__")
    __gt__ = __compare_by_reduce("__gt__")
    __le__ = __compare_by_reduce("__le__")
    __ge__ = __compare_by_reduce("__ge__")

    rtype = None

    @property
    def mol(self):
        """Get molecule.

        Returns:
            rdkit.Mol

        """
        return self._context.get_mol(self)

    @property
    def config(self):
        return self._context.config

    @property
    def coord(self):
        """Get 3D coordinate.

        Returns:
            numpy.array[3, N]: coordinate matrix

        """
        if not self.require_3D:
            self.fail(AttributeError("use 3D coordinate in 2D descriptor"))

        return self._context.get_coord(self)

    def get_3D_mol(self):
        mol = Chem.Mol(self.mol)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, xyz in enumerate(self.coord):
            conf.SetAtomPosition(i, xyz)

        mol.AddConformer(conf)
        return mol

    def fail(self, exception):
        """Raise known exception and return missing value.

        Raises:
            MissingValueException

        """
        raise MissingValueException(exception)

    @contextmanager
    def rethrow_zerodiv(self):
        """[contextmanager] treat zero div as known exception."""
        with np.errstate(divide="raise", invalid="raise"):
            try:
                yield
            except (FloatingPointError, ZeroDivisionError) as e:
                self.fail(ZeroDivisionError(*e.args))

    @contextmanager
    def rethrow_na(self, exception):
        """[contextmanager] treat any exceptions as known exception."""
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
                other = ConstDescriptor(other)

            return BinaryOperatingDescriptor(
                name.format(self, other), operator, self, other
            )

        return binary

    __add__ = _binary_common("({}+{})", "+")
    __sub__ = _binary_common("({}-{})", "-")
    __mul__ = _binary_common("({}*{})", "*")
    __truediv__ = _binary_common("({}/{})", "/")
    __floordiv__ = _binary_common("({}//{})", "//")
    __mod__ = _binary_common("({}%{})", "%")
    __pow__ = _binary_common("({}**{})", "**")

    __neg__ = _unary_common("-{}", "-")
    __pos__ = _unary_common("+{}", "+")
    __abs__ = _unary_common("|{}|", "abs")
    __trunc__ = _unary_common("trunc({})", "trunc")

    if six.PY3:
        __ceil__ = _unary_common("ceil({})", "ceil")
        __floor__ = _unary_common("floor({})", "floor")


def is_descriptor_class(desc, include_abstract=False):
    r"""Check calculatable descriptor class or not.

    Returns:
        bool

    """
    return (
        isinstance(desc, type)
        and issubclass(desc, Descriptor)
        and (True if include_abstract else not inspect.isabstract(desc))
    )


class UnaryOperatingDescriptor(Descriptor):
    @classmethod
    def preset(cls, version):
        return cls()

    operators = {
        "+": operator.pos,
        "-": operator.neg,
        "abs": operator.abs,
        "trunc": np.trunc,
        "ceil": np.ceil,  # noqa: S001
        "floor": np.floor,
    }

    def parameters(self):
        return self._name, self._operator, self._value

    def __init__(self, name, operator, value):
        self._name = name
        self._operator = operator
        self._fn = self.operators[operator]
        self._value = value

    def _to_json(self):
        return (
            self.__class__.__name__,
            {
                "name": self._name,
                "operator": self._operator,
                "value": self._value.to_json(),
            },
        )

    def __str__(self):
        return self._name

    def dependencies(self):
        return {"value": self._value}

    def calculate(self, value):
        return self._fn(value)


class ConstDescriptor(Descriptor):
    @classmethod
    def preset(cls, version):
        return cls()

    def parameters(self):
        return (self._value,)

    def __init__(self, value):
        self._value = value

    def __str__(self):
        return str(self._value)

    def calculate(self):
        return self._value


class BinaryOperatingDescriptor(Descriptor):
    @classmethod
    def preset(cls, version):
        return cls()

    operators = {
        "+": operator.add,
        "-": operator.sub,
        "*": operator.mul,  # noqa: S001
        "/": operator.truediv,
        "//": operator.floordiv,
        "%": operator.mod,  # noqa: S001
        "**": operator.pow,
    }

    def _to_json(self):
        return (
            self.__class__.__name__,
            {
                "name": self._name,
                "operator": self._operator,
                "left": self._left.to_json(),  # noqa: S001
                "right": self._right.to_json(),
            },
        )

    def parameters(self):
        return self._name, self._operator, self._left, self._right

    def __init__(self, name, operator, left, right):
        self._name = name
        self._operator = operator
        self._fn = self.operators[operator]
        self._left = left
        self._right = right

    def __str__(self):
        return self._name

    def dependencies(self):
        return {"left": self._left, "right": self._right}

    def calculate(self, left, right):
        return self._fn(left, right)
