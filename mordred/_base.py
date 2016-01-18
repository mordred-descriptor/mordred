from rdkit import Chem
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
from six import with_metaclass
from abc import ABCMeta, abstractmethod
from types import ModuleType
import os
from importlib import import_module
import numpy as np

from inspect import getsourcelines
from sys import maxsize


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


class Descriptor(with_metaclass(ABCMeta, object)):
    '''
    abstruct base class of descriptors
    '''

    explicit_hydrogens = True
    gasteiger_charges = False
    kekulize = False
    require_connected = True

    descriptor_keys = ()

    def __reduce_ex__(self, version):
        return self.__class__,\
            tuple(map(lambda k: getattr(self, k), self.descriptor_keys))

    def __hash__(self):
        return hash(tuple(map(lambda k: getattr(self, k), self.descriptor_keys)))

    def __eq__(self, other):
        return\
            self.__class__ is other.__class__ and\
            all(getattr(self, k) == getattr(other, k) for k in self.descriptor_keys)

    def __lt__(self, other):
        sk = tuple([self.__class__] + [getattr(self, k) for k in self.descriptor_keys])
        ok = tuple([other.__class__] + [getattr(other, k) for k in other.descriptor_keys])
        return sk.__lt__(ok)

    def __ne__(self, other):
        return not self == other

    @classmethod
    def preset(cls):
        '''
        generate preset descriptor instances

        Returns:
            iterable: descriptors
        '''
        yield cls()

    @property
    def dependencies(self):
        return None

    @abstractmethod
    def calculate(self, mol):
        pass

    def __call__(self, mol):
        '''
        calculate single descriptor value

        Returns:
            scalar: descriptor result
        '''
        return next(Calculator(self)(mol))[1]


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
    '''
    descriptor calculator

    Parameters:
        descs: see `register`
    '''

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
        '''
        register descriptors

        Parameters:
            descs (module, descriptor class/instance, iterable): descriptors to register
        '''

        for desc in descs:
            if not hasattr(desc, '__iter__'):
                if isinstance(desc, type) and issubclass(desc, Descriptor):
                    for d in desc.preset():
                        self._register_one(d)

                elif isinstance(desc, ModuleType):
                    self.register(self.get_descriptors_from_module(desc))

                else:
                    self._register_one(desc)

            else:
                for d in desc:
                    self.register(d)

    @staticmethod
    def get_descriptors_from_module(mdl):
        descs = []

        for name in dir(mdl):
            if name[:1] == '_':
                continue

            desc = getattr(mdl, name)
            if issubclass(desc, Descriptor):
                descs.append(desc)

        def key_by_def(d):
            try:
                return getsourcelines(d)[1]
            except IOError:
                return maxsize

        descs.sort(key=key_by_def)
        return descs

    def _calculate(self, desc, cache, parent=None):
        if desc in cache:
            return cache[desc]

        if desc.require_connected and not self.molecule.is_connected:
            cache[desc] = np.nan
            return np.nan

        args = {
            name: self._calculate(dep, cache, parent or desc)
            if dep is not None else None
            for name, dep in (desc.dependencies or {}).items()
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
        '''
        calculate descriptors

        Parameters:
            mol(rdkit.Chem.Mol): molecular

        Returns:
            iterator of descriptor and value

        :rtype: iterator of (Descriptor, scalar)
        '''
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
        u'''
        calculate descriptors over mols

        Parameters:
            mols(iterable<rdkit.Chem.Mol>): moleculars
            processes(int or None): number of process. None is multiprocessing.cpu_count()

        :rtype: iterator((rdkit.Chem.Mol, [(Descriptor, scalar)]]))
        '''

        if processes == 1:
            return ((m, list(self(m))) for m in mols)
        else:
            return self._parallel(mols, processes)


def initializer(calc):
    global calculate
    calculate = calc


def worker(binary):
    return list(calculate(Chem.Mol(binary)))


def all_descriptors():
    r'''
    yield all descriptors

    Returns:
        iterator<module>: all modules
    '''

    base_dir = os.path.dirname(__file__)

    for name in os.listdir(base_dir):
        name, ext = os.path.splitext(name)
        if name[:1] == '_' or ext != '.py':
            continue

        yield import_module('..' + name, __name__)
