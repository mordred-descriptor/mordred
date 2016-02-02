examples
--------

as command
~~~~~~~~~~

all descriptors
^^^^^^^^^^^^^^^

.. code:: console

    usage: python -m mordred [-h] [-i PATH] [-f TYPE] [-o PATH] [-p N] [-q] [-s]
                             [-3]

    optional arguments:
      -h, --help            show this help message and exit
      -i PATH, --input PATH
                            input file or directory(default: stdin)
      -f TYPE, --from TYPE  input filetype(one of auto, smi, sdf, mol, default:
                            auto)
      -o PATH, --output PATH
                            output csv file(default: stdout)
      -p N, --processes N   number of processes to use(default: number of threads)
      -q, --quiet           hide progress bar
      -s, --stream          stream read
      -3, --with-3D         calculate 3D descriptor(require sdf or mol file)

descriptors in submodule
^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: console

    $ python -m mordred.TPSA -i tests/references/structures.smi
    name,TPSA(NO),TPSA
    Hexane,0.0,0.0
    Benzene,0.0,0.0
    Caffeine,61.82,61.82
    Cyanidin,112.45,112.45
    Lycopene,0.0,0.0
    Epicatechin,110.38,110.38
    Limonene,0.0,0.0
    Allicin,17.07,61.58
    Glutathione,158.82,197.62

as library
^^^^^^^^^^

.. code:: python

    from rdkit import Chem

    from mordred import Calculator, all_descriptors

    # create descriptor calculator with all descriptors
    calc = Calculator(all_descriptors(with_3D=False))

    # calculate and print descriptors
    for desc, value in calc(Chem.MolFromSmiles('c1ccccc1O')):
       print('{}\t{}'.format(desc, value))

