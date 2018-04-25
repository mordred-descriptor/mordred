import inspect

from nose.tools import ok_

from mordred import descriptors
from mordred._base.descriptor import is_descriptor_class


def has__slots__(cls):
    return "__slots__" in cls.__dict__


def test_slots():
    for mdl in descriptors.all:
        for desc in mdl.__dict__.values():
            if not is_descriptor_class(desc, include_abstract=True):
                continue

            yield (
                ok_,
                has__slots__(desc),
                "{}({}) class don't have __slots__".format(
                    desc.__name__, inspect.getfile(desc),
                ),
            )
