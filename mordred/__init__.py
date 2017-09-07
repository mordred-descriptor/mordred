r"""modred descriptor calculator."""

from ._base import (
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
    "get_descriptors_from_module",
    "is_missing",
    "Result",
)
