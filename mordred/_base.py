import os
from abc import ABCMeta, abstractmethod
from importlib import import_module
from inspect import getsourcelines, isabstract
from sys import maxsize
from types import ModuleType

import numpy as np

from rdkit import Chem
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges

from six import with_metaclass


class DescriptorException(Exception):
    def __init__(self, desc, e, parent):
        self.desc = desc
        self.e = e
        self.parent = parent

    def __str__(self):
        if self.parent is None:
            return '{}.{}: {}'.format(
                self.desc.__class__.__name__,
                self.desc,
                self.e
            )

        return '{}.{}({}): {}'.format(
            self.desc.__class__.__name__,
            self.desc,
            self.parent,
            self.e
        )


def pretty(a):
    p = getattr(a, 'name', None)
    return repr(a if p is None else p)


class Descriptor(with_metaclass(ABCMeta, object)):
    r"""abstruct base class of descriptors."""

    explicit_hydrogens = True
    gasteiger_charges = False
    kekulize = False
    require_connected = True

    def _get_descriptor_keys(self):
        return (
            getattr(self, 'descriptor_keys', None) or
            getattr(self, '__slots__', None) or
            ()
        )

    def _get_keys(self):
        def getter(k):
            try:
                return getattr(self, k)
            except AttributeError as e:
                raise DescriptorException(self, e, None)

        return (getter(k) for k in self._get_descriptor_keys())

    def __reduce_ex__(self, version):
        return self.__class__, tuple(self._get_keys())

    def __repr__(self):
        return '{}({})'.format(
            self.__class__.__name__,
            ', '.join(map(pretty, self._get_keys()))
        )

    def __hash__(self):
        return hash(tuple(self._get_keys()))

    def __eq__(self, other):
        return\
            self.__class__ is other.__class__ and\
            all(getattr(self, k) == getattr(other, k) for k in self._get_descriptor_keys())

    def __lt__(self, other):
        sk = self.__reduce_ex__(3)
        ok = other.__reduce_ex__(3)
        return sk.__lt__(ok)

    def __ne__(self, other):
        return not self == other

    @classmethod
    def preset(cls):
        r"""generate preset descriptor instances.

        :rtype: iterable
        """
        return ()

    def dependencies(self):
        r"""descriptor dependencies.

        :rtype: {str: (Descriptor or None)} or None
        """
        return None

    @abstractmethod
    def calculate(self, mol):
        r"""calculate descriptor value.

        (abstruct method)
        """
        pass

    def __call__(self, mol):
        r"""calculate single descriptor value.

        :returns: descriptor result
        :rtype: scalar
        """
        return next(Calculator(self)(mol))[1]

    @classmethod
    def is_descriptor(cls, desc):
        return (
            isinstance(desc, type) and
            issubclass(desc, cls) and
            not isabstract(desc)
        )


class Molecule(object):
    def __init__(self, orig):
        Chem.SanitizeMol(orig)
        self.orig = orig
        self.hydration_cache = dict()
        self.kekulize_cache = dict()
        self.gasteiger_cache = dict()
        self.is_connected = len(Chem.GetMolFrags(orig)) == 1

    def hydration(self, explicitH):
        if explicitH in self.hydration_cache:
            return self.hydration_cache[explicitH]

        mol = Chem.AddHs(self.orig) if explicitH else Chem.RemoveHs(self.orig)
        self.hydration_cache[explicitH] = mol
        return mol

    def kekulize(self, mol, explicitH):
        if explicitH in self.kekulize_cache:
            return self.kekulize_cache[explicitH]

        mol = Chem.Mol(mol)
        Chem.Kekulize(mol)
        self.kekulize_cache[explicitH] = mol
        return mol

    def gasteiger(self, mol, explicitH, kekulize):
        key = explicitH, kekulize
        if key in self.gasteiger_cache:
            return self.gasteiger_cache[key]

        ComputeGasteigerCharges(mol)
        self.gasteiger_cache[key] = mol
        return mol

    def get(self, explicitH, kekulize, gasteiger):
        mol = self.hydration(explicitH)
        if kekulize:
            mol = self.kekulize(mol, explicitH)
        if gasteiger:
            mol = self.gasteiger(mol, explicitH, kekulize)

        return mol


class Calculator(object):
    r"""descriptor calculator.

    :param descs: see `register` method
    """

    def __init__(self, *descs):
        self.descriptors = []
        self.explicitH = False
        self.gasteiger = False
        self.kekulize = False

        self.register(*descs)

    def __reduce_ex__(self, version):
        return self.__class__, tuple(self.descriptors)

    def _register_one(self, desc):
        if not isinstance(desc, Descriptor):
            raise ValueError('{!r} is not descriptor'.format(desc))

        self.descriptors.append(desc)

        if desc.explicit_hydrogens:
            self.explicitH = True

        if desc.gasteiger_charges:
            self.gasteiger = True

        if desc.kekulize:
            self.kekulize = True

    def register(self, *descs):
        r"""register descriptors.

        :param descs: descriptors to register
        :type descs: module, descriptor class/instance, iterable
        """
        for desc in descs:
            if not hasattr(desc, '__iter__'):
                if Descriptor.is_descriptor(desc):
                    for d in desc.preset():
                        self._register_one(d)

                elif isinstance(desc, ModuleType):
                    self.register(get_descriptors_from_module(desc))

                else:
                    self._register_one(desc)

            else:
                for d in desc:
                    self.register(d)

    def _calculate(self, desc, cache, parent=None):
        if desc in cache:
            return cache[desc]

        if desc.require_connected and not self.molecule.is_connected:
            cache[desc] = np.nan
            return np.nan

        args = {
            name: self._calculate(dep, cache, parent or desc)
            if dep is not None else None
            for name, dep in (desc.dependencies() or {}).items()
        }

        mol = self.molecule.get(
            explicitH=desc.explicit_hydrogens,
            gasteiger=desc.gasteiger_charges,
            kekulize=desc.kekulize,
        )

        try:
            r = desc.calculate(mol, **args)
        except Exception as e:
            raise DescriptorException(desc, e, parent)

        cache[desc] = r
        return r

    def __call__(self, mol):
        r"""calculate descriptors.

        Parameters:
            mol(rdkit.Chem.Mol): molecular

        Returns:
            iterator of descriptor and value

        :rtype: iterator of (Descriptor, scalar)
        """
        cache = {}
        self.molecule = Molecule(mol)

        return (
            (desc, self._calculate(desc, cache))
            for desc in self.descriptors
        )

    def _parallel(self, mols, processes=None):
        from multiprocessing import Pool

        try:
            pool = Pool(
                processes,
                initializer=initializer,
                initargs=(self,),
            )

            for m, result in [(m, pool.apply_async(worker, (m.ToBinary(),))) for m in mols]:
                yield m, result.get()

        finally:
            pool.terminate()

    def map(self, mols, processes=None):
        r"""calculate descriptors over mols.

        :param mols: moleculars
        :type mols: iterable(rdkit.Chem.Mol)

        :param processes: number of process. None is multiprocessing.cpu_count()
        :type processes: int or None

        :rtype: iterator((rdkit.Chem.Mol, [(Descriptor, scalar)]]))
        """
        if processes == 1:
            return ((m, list(self(m))) for m in mols)
        else:
            return self._parallel(mols, processes)


calculate = None


def initializer(calc):
    global calculate
    calculate = calc


def worker(binary):
    return list(calculate(Chem.Mol(binary)))


def all_descriptors():
    r"""yield all descriptors.

    :returns: all modules
    :rtype: iterator(module)
    """
    base_dir = os.path.dirname(__file__)

    for name in os.listdir(base_dir):
        name, ext = os.path.splitext(name)
        if name[:1] == '_' or ext != '.py':
            continue

        yield import_module('..' + name, __name__)


def get_descriptors_from_module(mdl):
    r"""get descriptors from module.

    :type mdl: module

    :rtype: [Descriptor]
    """
    descs = []

    for name in dir(mdl):
        if name[:1] == '_':
            continue

        desc = getattr(mdl, name)
        if Descriptor.is_descriptor(desc):
            descs.append(desc)

    def key_by_def(d):
        try:
            return getsourcelines(d)[1]
        except IOError:
            return maxsize

    descs.sort(key=key_by_def)
    return descs


def parse_enum(enum, v):
    if isinstance(v, enum):
        return v
    else:
        return enum[v]
