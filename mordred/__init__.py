r"""modred descriptor calculator."""

from ._version import __version__

from ._base import (
    all_descriptors,
    Calculator,
    Descriptor,
    get_descriptors_from_module,
    is_missing,
)

__all__ = (
    '__version__',
    'Descriptor',
    'Calculator',
    'all_descriptors',
    'get_descriptors_from_module',
    'is_missing',
)
