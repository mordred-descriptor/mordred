import load_path
from mordred import Descriptor, descriptors, get_descriptors_from_module

load_path.nop()

prelude = '''
Descriptor List
===============
preset descriptors

.. code:: python

    calc = Calculator(descriptors)

.. csv-table:: Descriptor list
    :header: "#", "module", "name", "constructor", "dim"
    :widths: 10, 20, 20, 40, 10

'''[1:]


def main(out):
    out.write(prelude)

    i = 0

    for mdl in descriptors.all:
        mdl_name = '.'.join(mdl.__name__.split('.'))
        mdl_ppr = ':py:mod:`~{}`'.format(mdl_name)
        first = True

        for Desc in get_descriptors_from_module(mdl):

            for desc in Desc.preset():
                i += 1

                if not first:
                    mdl_ppr = ''

                cnst = desc.__class__.__name__
                args = ', '.join(Descriptor._pretty(p) for p in desc.parameters())

                cnst = ':py:class:`~{}.{}` ({})'.format(mdl_name, cnst, args)

                dim = '3D' if desc.require_3D else '2D'

                out.write('    {}, {}, {}, "{}", {}\n'.format(i, mdl_ppr, desc, cnst, dim))

                first = False


if __name__ == '__main__':
    import sys
    main(sys.stdout)
