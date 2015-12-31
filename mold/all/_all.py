import os
from glob import glob
from importlib import import_module


__all__ = ['descriptors']


desc_pathes = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    '*', '_descriptor.py'
)


descriptors = []
for path in glob(desc_pathes):
    desc_name = os.path.basename(os.path.dirname(path))
    if desc_name[:1] == '_':
        continue

    mdl = import_module('...' + desc_name, __name__)
    descriptors.append(mdl)

    globals()[desc_name] = mdl
    __all__.append(desc_name)
