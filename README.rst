mordred
=======
molecular descriptor calculator.

.. image:: https://travis-ci.org/mordred-descriptor/mordred.svg?branch=master
    :target: https://travis-ci.org/mordred-descriptor/mordred

.. image:: https://ci.appveyor.com/api/projects/status/iwk6268d27jusvni/branch/master?svg=true
    :target: https://ci.appveyor.com/project/philopon/mordred/branch/master
    
.. image:: https://coveralls.io/repos/github/mordred-descriptor/mordred/badge.svg?branch=master
    :target: https://coveralls.io/github/mordred-descriptor/mordred?branch=master 

.. image:: https://codeclimate.com/github/mordred-descriptor/mordred/badges/gpa.svg
   :target: https://codeclimate.com/github/mordred-descriptor/mordred
   :alt: Code Climate

.. image:: https://anaconda.org/mordred-descriptor/mordred/badges/version.svg
    :target: https://anaconda.org/mordred-descriptor/mordred

Installation
------------

conda(recommended)
~~~~~~~~~~~~~~~~~~
#. install conda

       -  `miniconda <http://conda.pydata.org/miniconda.html>`__
       -  `anaconda <https://www.continuum.io/why-anaconda>`__

#. install mordred

       stable

       .. code:: console

           $ conda install -c rdkit -c mordred-descriptor mordred

       beta

       .. code:: console

           $ conda install -c rdkit -c mordred-descriptor/channel/dev mordred

pip
~~~

#. install `rdkit <http://www.rdkit.org/>`__ python package
#. install mordred

       .. code:: console

           $ pip install git+https://github.com/mordred-descriptor/mordred

examples
--------

as command
~~~~~~~~~~

.. code:: console

    usage: python -m mordred [-h] [-f TYPE] [-o PATH] [-p N] [-q] [-s] [-3] INPUT

    positional arguments:
      INPUT                 input file or directory(default: stdin)

    optional arguments:
      -h, --help            show this help message and exit
      -f TYPE, --from TYPE  input filetype(one of auto, smi, sdf, mol, default:
                            auto)
      -o PATH, --output PATH
                            output csv file(default: stdout)
      -p N, --processes N   number of processes to use(default: number of threads)
      -q, --quiet           hide progress bar
      -s, --stream          stream read
      -3, --with-3D         calculate 3D descriptor(require sdf or mol file)

as library
^^^^^^^^^^

.. code:: python

    from rdkit import Chem

    from mordred import Calculator, all_descriptors

    # create descriptor calculator with all descriptors
    calc = Calculator(all_descriptors(), exclude3D=True)

    mol = Chem.MolFromSmiles('c1ccccc1')

    mols = map(Chem.MolFromSmiles, [
        'c1ccccc1Cl',
        'c1ccccc1O',
        'c1ccccc1N'
    ])

    # calculate single molecule
    for desc, value in zip(calc.descriptors, calc(mol)):
        print('{}\t{}'.format(desc, value))

    # calculate multiple molecule
    for mol, values in calc.map(mols, processes=1):
        print(Chem.MolToSmiles(mol))
        for desc, value in zip(calc.descriptors, values):
            print('{}\t{}'.format(desc, value))

Documentation
-------------

-  `v0.2.0 <http://mordred-descriptor.github.io/documentation/v0.2.0>`__
-  `v0.1.0 <http://mordred-descriptor.github.io/documentation/v0.1.0>`__

