r"""modred descriptor calculator."""

from ._base import (
    Result,
    Calculator,
    Descriptor,
    is_missing,
    get_descriptors_in_module,
    get_descriptors_from_module,
)
from ._version import __version__

__all__ = (
    "__version__",
    "Descriptor",
    "Calculator",
    "get_descriptors_from_module",
    "get_descriptors_in_module",
    "is_missing",
    "Result",
)
