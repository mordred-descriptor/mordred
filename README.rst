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

.. image:: https://img.shields.io/pypi/v/mordred.svg
    :target: https://pypi.python.org/pypi/mordred

.. image:: https://img.shields.io/badge/doi-10.1186%2Fs13321--018--0258--y-blue.svg
   :target: https://doi.org/10.1186/s13321-018-0258-y


Installation
------------

conda(recommended)
~~~~~~~~~~~~~~~~~~
#. install conda

       -  `miniconda <http://conda.pydata.org/miniconda.html>`__
       -  `anaconda <https://www.continuum.io/why-anaconda>`__

#. install mordred

       .. code:: console

           $ conda install -c rdkit -c mordred-descriptor mordred

pip
~~~

#. install `rdkit <http://www.rdkit.org/>`__ python package
#. install mordred

       .. code:: console

           $ pip install mordred

examples
--------

as command
~~~~~~~~~~

calculate all descriptors

.. code:: console

    $ python -m mordred example.smi
    name,ECIndex,WPath,WPol,Zagreb1, (snip)
    benzene,36,27,3,24.0, (snip)
    chrolobenzene,45,42,5,30.0, (snip)


save to file (display progress bar)

.. code:: console

    $ python -m mordred example.smi -o example.csv
    50%|███████████████████████████████████████▌                                       | 1/2 [00:00<00:00,  7.66it/s]


stream read (low memory, no number of molecules information)

.. code:: console

    $ python -m mordred example.smi -s -o example.csv
    0it [00:00, ?it/s]

only ABCIndex

.. code:: console

    $ python -m mordred example.smi -d ABCIndex
    name,ABC,ABCGG
    benzene,4.242640687119286,3.9999999999999996
    chlorobenzene,5.059137268047012,4.785854275382693

ABCIndex and AcidBase

.. code:: console

    $ python -m mordred example.smi -d ABCIndex -d AcidBase
    name,ABC,ABCGG,nAcid,nBase
    benzene,4.242640687119286,3.9999999999999996,0,0
    chlorobenzene,5.059137268047012,4.785854275382693,0,0

multiple input

.. code:: console

    $ python -m mordred example.smi example2.smi -d ABCIndex
    name,ABC,ABCGG
    benzene,4.242640687119286,3.9999999999999996
    chlorobenzene,5.059137268047012,4.785854275382693
    pentane,2.8284271247461903,3.1462643699419726

show help

.. code:: console

    $ python -m mordred --help
    usage: python -m mordred [-h] [--version] [-t {auto,sdf,mol,smi}] [-o OUTPUT]
                             [-p PROCESSES] [-q] [-s] [-d DESC] [-3] [-v]
                             INPUT [INPUT ...]

    positional arguments:
      INPUT

    optional arguments:
      -h, --help            show this help message and exit
      --version             input molecular file
      -t {auto,sdf,mol,smi}, --type {auto,sdf,mol,smi}
                            input filetype (default: auto)
      -o OUTPUT, --output OUTPUT
                            output file path (default: stdout)
      -p PROCESSES, --processes PROCESSES
                            number of processes (default: number of logical
                            processors)
      -q, --quiet           hide progress bar
      -s, --stream          stream read
      -d DESC, --descriptor DESC
                            descriptors to calculate (default: all)
      -3, --3D              use 3D descriptors (require sdf or mol file)
      -v, --verbosity       verbosity

    descriptors: ABCIndex AcidBase AdjacencyMatrix Aromatic AtomCount
    Autocorrelation BalabanJ BaryszMatrix BCUT BertzCT BondCount CarbonTypes Chi
    Constitutional CPSA DetourMatrix DistanceMatrix EccentricConnectivityIndex
    EState ExtendedTopochemicalAtom FragmentComplexity Framework GeometricalIndex
    GravitationalIndex HydrogenBond InformationContent KappaShapeIndex Lipinski
    McGowanVolume MoeType MolecularDistanceEdge MolecularId MomentOfInertia MoRSE
    PathCount Polarizability RingCount RotatableBond SLogP TopologicalCharge
    TopologicalIndex TopoPSA VdwVolumeABC VertexAdjacencyInformation WalkCount
    Weight WienerIndex ZagrebIndex

as library
^^^^^^^^^^

.. code:: python

    >>> from rdkit import Chem
    >>> from mordred import Calculator, descriptors

    # create descriptor calculator with all descriptors
    >>> calc = Calculator(descriptors, ignore_3D=True)

    >>> len(calc.descriptors)
    1612

    # calculate single molecule
    >>> mol = Chem.MolFromSmiles('c1ccccc1')
    >>> calc(mol)[:3]
    [4.242640687119286, 3.9999999999999996, 0]

    # calculate multiple molecule
    >>> mols = [Chem.MolFromSmiles(smi) for smi in ['c1ccccc1Cl', 'c1ccccc1O', 'c1ccccc1N']]

    # as pandas
    >>> df = calc.pandas(mols)
    >>> df['SLogP']
    0    2.3400
    1    1.3922
    2    1.2688
    Name: SLogP, dtype: float64

see `examples <examples>`_

Citation
--------
Moriwaki H, Tian Y-S, Kawashita N, Takagi T (2018) Mordred: a molecular descriptor calculator. Journal of Cheminformatics 10:4 . doi: `10.1186/s13321-018-0258-y <https://doi.org/10.1186/s13321-018-0258-y>`__

Documentation
-------------

-  `master <http://mordred-descriptor.github.io/documentation/master>`__
-  `develop <http://mordred-descriptor.github.io/documentation/develop>`__

-  `v1.0.0 <http://mordred-descriptor.github.io/documentation/v1.0.0>`__
