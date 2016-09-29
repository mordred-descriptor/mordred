import os
from importlib import import_module
from ..error import MissingValueBase

from .descriptor import Descriptor
from .calculator import Calculator, get_descriptors_from_module
from .parallel import parallel


__all__ = (
    'all_descriptors',
    'Descriptor',
    'Calculator',
    'get_descriptors_from_module',
)


def all_descriptors():
    r"""yield all descriptor modules.

    :returns: all modules
    :rtype: :py:class:`Iterator` (:py:class:`Descriptor`)
    """
    base_dir = os.path.dirname(os.path.dirname(__file__))

    for name in os.listdir(base_dir):
        name, ext = os.path.splitext(name)
        if name[:1] == '_' or ext != '.py':
            continue

        yield import_module('..' + name, __package__)


def Descriptor__call__(self, mol, id=-1):
    r"""calculate single descriptor value.
    :type id: int
    :param id: conformer id

    :returns: descriptor result
    :rtype: scalar
    """
    v = Calculator(self)(mol, id)[0]
    if isinstance(v, MissingValueBase):
        raise v.error

    return v


Descriptor.__call__ = Descriptor__call__
Calculator._parallel = parallel
