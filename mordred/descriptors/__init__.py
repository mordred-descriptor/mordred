"""All descriptor modules are loaded.

.. code:: python

    >>> from mordred import Calculator, descriptors

    >>> descriptors.ABCIndex  # ABCIndex module
    <module 'mordred.ABCIndex' from '/path/to/ABCIndex.py'>

    >>> descriptors.all  # all descriptor modules
    (<module 'mordred.ABCIndex' from '/path/to/ABCIndex.py'>,
     <module 'mordred.AcidBase' from '/path/to/AcidBase.py'>,
     ...
    )

    >>> calc = Calculator(descriptors) # all descriptors
"""


def _import_all_descriptors():
    import os
    from importlib import import_module
    from .._base.descriptor import is_descriptor_class

    names = []
    values = []
    base_dir = os.path.dirname(os.path.dirname(__file__))

    for name in os.listdir(base_dir):
        name, ext = os.path.splitext(name)
        if name[:1] == "_" or ext != ".py":
            continue

        mdl = import_module(".." + name, __package__)

        if any(v for v in mdl.__dict__.values() if is_descriptor_class(v)):
            names.append(name)
            values.append(mdl)
            globals()[name] = mdl

        globals()["__all__"] = tuple(names)
        globals()["all"] = tuple(values)


_import_all_descriptors()
