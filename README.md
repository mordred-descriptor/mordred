# mordred-community
Community maintained version of the [mordred molecular descriptor calculator](https://github.com/mordred-descriptor/mordred) which is no longer maintained.

We are accepting pull requests and looking for maintainers!
Reach out on the Issues page if you are interested in helping out!

## Installation
`mordredcommunity` is currently available primarily on PyPI (conda package (_hopefully_) coming soon).

It supports Python 3.7 to 3.11 on all platforms.

For a basic installation, do `pip install mordredcommunity`.

To add support for `pandas` and progress bars, do `pip install mordredcommunity[full]`

Usage of `mordredcommunity` is the same as `mordred`.

## Number of Descriptors

```python
>>> from mordred import Calculator, descriptors
>>> n_all = len(Calculator(descriptors, ignore_3D=False).descriptors)
>>> n_2D = len(Calculator(descriptors, ignore_3D=True).descriptors)
>>> print("2D:    {:5}\n3D:    {:5}\n------------\ntotal: {:5}".format(n_2D, n_all - n_2D, n_all))
2D:     1613
3D:      213
------------
total:  1826
```

# Examples

## calculate all descriptors

```python
$ python -m mordred example.smi
name,ECIndex,WPath,WPol,Zagreb1, (snip)
benzene,36,27,3,24.0, (snip)
chrolobenzene,45,42,5,30.0, (snip)
```


## save to file (display progress bar)

```python
$ python -m mordred example.smi -o example.csv
50%|███████████████████████████████████████▌                                       | 1/2 [00:00<00:00,  7.66it/s]
```

## stream read (low memory, no number of molecules information)

```python
$ python -m mordred example.smi -s -o example.csv
0it [00:00, ?it/s]
```

## only ABCIndex

```python
$ python -m mordred example.smi -d ABCIndex
name,ABC,ABCGG
benzene,4.242640687119286,3.9999999999999996
chlorobenzene,5.059137268047012,4.785854275382693
```

## ABCIndex and AcidBase

```python
$ python -m mordred example.smi -d ABCIndex -d AcidBase
name,ABC,ABCGG,nAcid,nBase
benzene,4.242640687119286,3.9999999999999996,0,0
chlorobenzene,5.059137268047012,4.785854275382693,0,0
```

## multiple input

```python
$ python -m mordred example.smi example2.smi -d ABCIndex
name,ABC,ABCGG
benzene,4.242640687119286,3.9999999999999996
chlorobenzene,5.059137268047012,4.785854275382693
pentane,2.8284271247461903,3.1462643699419726
```

## show help

```python
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
```

## as library

```python
>>> from rdkit import Chem
>>> from mordred import Calculator, descriptors

# create descriptor calculator with all descriptors
>>> calc = Calculator(descriptors, ignore_3D=True)

>>> len(calc.descriptors)
1613

>>> len(Calculator(descriptors, ignore_3D=True, version="1.0.0"))
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
```

## Additional Examples
see [examples](https://github.com/JacksonBurns/mordred-community/tree/main/examples) on GitHub.

# Citation
Please cite the original publication describing `mordred` _and_ this repository specifically.
Moriwaki H, Tian Y-S, Kawashita N, Takagi T (2018) Mordred: a molecular descriptor calculator. Journal of Cheminformatics 10:4 . doi: [10.1186/s13321-018-0258-y](https://doi.org/10.1186/s13321-018-0258-y)
