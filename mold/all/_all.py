import os
from glob import glob
from importlib import import_module


__all__ = ['descriptors']


descriptors = []

for path in glob(os.path.join(os.path.dirname(__file__), '..', '*', '_descriptor.py')):
    desc_name = os.path.basename(os.path.dirname(path))
    if desc_name[:1] == '_':
        continue

    mdl = import_module('...' + desc_name, __name__)

    desc_added = False
    if hasattr(mdl, '_descriptors'):
        descriptors += mdl._descriptors
        desc_added = True

    for name in dir(mdl):
        if name[:1] != '_':
            desc = getattr(mdl, name)
            globals()[name] = desc
            __all__.append(name)

            if not desc_added:
                descriptors.append(desc)
