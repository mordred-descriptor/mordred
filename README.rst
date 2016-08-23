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

       .. code:: console

           $ conda install -c rdkit -c mordred-descriptor mordred

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
    Autocorrelation Autocorrelation3D BalabanJ BaryszMatrix BCUT BertzCT BondCount
    CarbonTypes Chi Constitutional CPSA DetourMatrix DistanceMatrix
    EccentricConnectivityIndex EState ExtendedTopochemicalAtom FragmentComplexity
    Framework GeometricalIndex GravitationalIndex HydrogenBond InformationContent
    KappaShapeIndex Lipinski McGowanVolume MoeType MolecularDistanceEdge
    MolecularId MomentOfInertia MoRSE PathCount Polarizability RingCount
    RotatableBond SLogP TopologicalCharge TopologicalIndex TopoPSA VdwVolumeABC

as library
^^^^^^^^^^

.. code:: python

    >>> from rdkit import Chem
    >>> from mordred import Calculator, all_descriptors

    # create descriptor calculator with all descriptors
    >>> calc = Calculator(all_descriptors(), exclude3D=True)

    >>> len(calc.descriptors)
    1824

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

Documentation
-------------

-  `v0.2.0 <http://mordred-descriptor.github.io/documentation/v0.2.0>`__
-  `v0.1.0 <http://mordred-descriptor.github.io/documentation/v0.1.0>`__

