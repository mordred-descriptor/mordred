"""Mordred base package."""

from .util import is_missing
from ..error import MissingValueBase
from .result import Result
from .parallel import parallel
from .calculator import (
    Calculator,
    get_descriptors_in_module,
)
from .descriptor import (
    Descriptor,
    ConstDescriptor,
    UnaryOperatingDescriptor,
    BinaryOperatingDescriptor,
)

__all__ = (
    "Descriptor",
    "Calculator",
    "get_descriptors_in_module",
    "is_missing",
    "Result",
)


def _Descriptor__call__(self, mol, id=-1):
    r"""Calculate single descriptor value.

    :type id: int
    :param id: conformer id

    :returns: descriptor result
    :rtype: scalar
    """
    v = Calculator(self)(mol, id)[0]
    if isinstance(v, MissingValueBase):
        raise v.error

    return v


def _from_json(obj, descs):
    name = obj.get("name")
    args = obj.get("args") or {}
    if name is None:
        raise ValueError("invalid json: {}".format(obj))

    if name == UnaryOperatingDescriptor.__name__:
        return UnaryOperatingDescriptor(
            args["name"], args["operator"], _from_json(args["value"])
        )

    elif name == BinaryOperatingDescriptor.__name__:
        return BinaryOperatingDescriptor(
            args["name"],
            args["operator"],
            _from_json(args["left"]),
            _from_json(args["right"]),
        )

    cls = descs.get(name)
    if cls is None:
        raise ValueError("unknown class: {}".format(name))

    instance = cls(**(obj.get("args") or {}))
    return instance


@classmethod
def _Descriptor_from_json(self, obj):
    """Create Descriptor instance from json dict.

    Parameters:
        obj(dict): descriptor dict

    Returns:
        Descriptor: descriptor

    """
    descs = getattr(self, "_all_descriptors", None)

    if descs is None:
        from mordred import descriptors

        descs = {cls.__name__: cls for cls in get_descriptors_in_module(descriptors)}
        descs[ConstDescriptor.__name__] = ConstDescriptor
        self._all_descriptors = descs

    return _from_json(obj, descs)


Descriptor.__call__ = _Descriptor__call__
Descriptor.from_json = _Descriptor_from_json
Calculator._parallel = parallel
