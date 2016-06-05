r"""modred descriptor calculator."""

from ._version import __version__

from ._base import (
    all_descriptors,
    Calculator, Descriptor, Error,
    get_descriptors_from_module,
)

__all__ = (
    '__version__',
    'Descriptor', 'Error',
    'Calculator',
    'all_descriptors',
    'get_descriptors_from_module',
)
