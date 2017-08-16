import os
import warnings
from importlib import import_module

from ..error import MissingValueBase


def all_descriptors():
    r"""**[deprecated]** use mordred.descriptors module instead.

    yield all descriptor modules.

    :returns: all modules
    :rtype: :py:class:`Iterator` (:py:class:`Descriptor`)
    """
    warnings.warn(
        "all_descriptors() is deprecated, use mordred.descriptors module instead",
        DeprecationWarning,
        stacklevel=2,
    )
    base_dir = os.path.dirname(os.path.dirname(__file__))

    for name in os.listdir(base_dir):
        name, ext = os.path.splitext(name)
        if name[:1] == "_" or ext != ".py" or name == "descriptors":
            continue

        yield import_module(".." + name, __package__)


def is_missing(v):
    """Check argument is either MissingValue or not.

    Parameters:
        v(any): value

    Returns:
        bool

    """
    return isinstance(v, MissingValueBase)
