mordred [![Build Status](https://travis-ci.org/mordred-descriptor/mordred.svg?branch=master)](https://travis-ci.org/mordred-descriptor/mordred) [![Coverage Status](https://coveralls.io/repos/mordred-descriptor/mordred/badge.svg?branch=master&service=github)](https://coveralls.io/github/mordred-descriptor/mordred?branch=master) [![Code Climate](https://codeclimate.com/github/mordred-descriptor/mordred/badges/gpa.svg)](https://codeclimate.com/github/mordred-descriptor/mordred) [![Anaconda-Server Badge](https://anaconda.org/mordred-descriptor/mordred/badges/version.svg)](https://anaconda.org/mordred-descriptor/mordred)
==
molecular descriptor calculator

Installation
--
### conda(recommended)

1. install conda

    * [miniconda](http://conda.pydata.org/miniconda.html)
    * [anaconda](https://www.continuum.io/why-anaconda)

2. install mordred

    #### stable
    
    ```
    $ conda install -c rdkit -c mordred-descriptor mordred
    ```
    
    #### beta
    
    ```
    $ conda install -c rdkit -c mordred-descriptor/channel/dev mordred
    ```

### pip

1. install [rdkit](http://www.rdkit.org/) python package

2. install mordred

   ```
   $ pip install git+https://github.com/mordred-descriptor/mordred
   ```

Usage
--

### as command

#### all descriptors

```console
$ python -m mordred --help
usage: mordred [-h] [-i [PATH]] [-f TYPE] [-o [PATH]] [-p N]

optional arguments:
  -h, --help            show this help message and exit
  -i [PATH], --input [PATH]
                        input file(default: stdin)
  -f TYPE, --from TYPE  input filetype(one of auto, smi, sdf, default: auto)
  -o [PATH], --output [PATH]
                        output csv file(default: stdout)
  -p N, --processes N   number of processes to use(default: number of threads)
```

#### descriptors in submodule

```console
$ python -m mordred.TPSA -i tests/data/structures.smi
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
```

### as library

```.py
from rdkit import Chem

from mordred import Calculator, all_descriptors

# create descriptor calculator with all descriptors
calc = Calculator(all_descriptors())

# calculate and print descriptors
for desc, value in calc(Chem.MolFromSmiles('c1ccccc1O')):
   print('{}\t{}'.format(desc, value))
```

Documentation
--
* [stable](http://mordred-descriptor.github.io/documentation/release)
* [beta](http://mordred-descriptor.github.io/documentation/master)
