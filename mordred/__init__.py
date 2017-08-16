r"""modred descriptor calculator."""

from ._base import (
    all_descriptors,
    Calculator,
    Descriptor,
    get_descriptors_from_module,
    is_missing,
    Result,
)

from ._version import __version__

__all__ = (
    "__version__",
    "Descriptor",
    "Calculator",
    "all_descriptors",
    "get_descriptors_from_module",
    "is_missing",
    "Result",
)
