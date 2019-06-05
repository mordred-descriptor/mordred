import sys
from distutils.version import StrictVersion

import load_path
from mordred import Descriptor, __version__, descriptors, get_descriptors_in_module

load_path.nop()

prelude = """
Descriptor List
===============
preset descriptors

.. code:: python

    calc = Calculator(descriptors)

.. csv-table:: Descriptor list
    :header: "#", "module", "name", "constructor", "dim", "description"
    :widths: 1, 2, 2, 4, 1, 10

"""[1:]


class DescriptorInfo(object):
    def __init__(self, d):
        self.raw = d

    @property
    def module(self):
        return self.raw.__module__

    @property
    def constructor(self):
        return self.raw.__class__.__name__

    @property
    def parameters(self):
        return [Descriptor._pretty(p) for p in self.raw.parameters()]

    @property
    def dimention(self):
        return "3D" if self.raw.require_3D else "2D"

    @property
    def description(self):
        return self.raw.description()

    def to_rst(self, hide_module=False):
        mdl = "" if hide_module else ":py:mod:`~{}`".format(self.module)
        desc = self.raw
        cnst = ":py:class:`~{}.{}` ({})".format(
            self.module, self.constructor, ", ".join(self.parameters))
        info = self.description
        dim = self.dimention

        return '{}, {}, "{}", {}, "{}"'.format(mdl, desc, cnst, dim, info)


def get_all_descriptors():
    v = StrictVersion(__version__)
    for mdl in descriptors.all:
        ds = []
        for Desc in get_descriptors_in_module(mdl, submodule=False):
            for desc in Desc.preset(v):
                ds.append(desc)

        yield ds


def main(out):
    out.write(prelude)

    i = 0

    for descs in get_all_descriptors():
        first = True
        for desc in descs:
            i += 1
            out.write("    {}, {}\n".format(i, DescriptorInfo(desc).to_rst(not first)))
            first = False


if __name__ == "__main__":
    main(sys.stdout)
