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

example
--

```.py
from rdkit import Chem

from mordred import Calculator
import mordred.all

# create descriptor calculator with all descriptors
calc = Calculator(mordred.all.descriptors)

# calculate and print descriptors
for name, value in calc(Chem.MolFromSmiles('c1ccccc1O')):
   print(name, value)
```

documentation
--
* [stable](http://mordred-descriptor.github.io/documentation/release)
* [beta](http://mordred-descriptor.github.io/documentation/master)
