import conf
from mordred import all_descriptors, get_descriptors_from_module

prelude = '''
Descriptor List
===============
preset descriptors

.. code:: python

    calc = Calculator(all_descriptors())

.. csv-table:: Descriptor list
    :header: "#", "module", "name", "constructor", "dim"
    :widths: 10, 20, 20, 40, 10

'''[1:]


def main(out):
    out.write(prelude)
    
    i = 0

    for mdl in all_descriptors():
        for Desc in get_descriptors_from_module(mdl):
            for desc in Desc.preset():
                i += 1

                if mdl:
                    mdl_name = '.'.join(mdl.__name__.split('.')[1:])
                    mdl_ppr = ':py:mod:`~mordred.{}`'.format(mdl_name)
                    mdl = None
                else:
                    mdl_ppr = ''

                try:
                    cnst, args = repr(desc).split('(')
                except ValueError:
                    raise ValueError(repr(desc))

                cnst = ':py:class:`~mordred.{}.{}` ({}'.format(mdl_name, cnst, args)

                dim = '3D' if desc.require_3D else '2D'

                out.write('    {}, {}, {}, "{}", {}\n'.format(i, mdl_ppr, desc, cnst, dim))


if __name__ == '__main__':
    import sys
    main(sys.stdout)
