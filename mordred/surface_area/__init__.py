r"""surface area subpackage.

.. code:: console

    $ python -m mordred.surface_area -h
    usage: python -m mordred.surface_area [-h] [-s R] [-l L] FILE

    positional arguments:
      FILE                  input sdf/mol file

    optional arguments:
      -h, --help            show this help message and exit
      -s R, --solvent-radius R
                            solvent radius (default: 1.4)
      -l L, --mesh-level L  mesh level (default: 5)
"""

from ._sasa import SurfaceArea

__all__ = ("SurfaceArea",)
