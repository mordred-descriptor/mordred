import inspect

from mordred import descriptors, get_descriptors_in_module
from mordred._base.descriptor import is_descriptor_class


def test_attributes():
    for mdl in descriptors.all:
        for desc in get_descriptors_in_module(mdl, submodule=False):
            for cls in desc.__mro__:
                if cls == object:
                    continue

                if is_descriptor_class(desc, include_abstract=True):
                    assert (
                        "__slots__" in cls.__dict__
                    ), "{}({}) class don't have __slots__".format(
                        cls.__name__, inspect.getfile(desc)
                    )

            assert "since" in desc.__dict__, "{}({}) class don't have since".format(
                desc.__name__, inspect.getfile(desc)
            )
