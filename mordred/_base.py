from rdkit import Chem
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
from six import with_metaclass
from abc import ABCMeta, abstractmethod
from types import ModuleType
import os
from importlib import import_module

from inspect import getsourcelines
from sys import maxsize

try:
    from inspect import getfullargspec as getargspec
except ImportError:
    from inspect import getargspec


class Key(object):
    __slots__ = 'cls', 'args'

    def __init__(self, cls, *args):
        assert isinstance(cls, type)
        self.cls = cls
        self.args = args

    def __hash__(self):
        return hash((self.cls, self.args))

    def __eq__(self, other):
        return self.cls is other.cls and self.args == other.args

    def __ne__(self, other):
        return self.cls is not other.cls or self.args != other.args

    def __str__(self):
        return self.cls.__name__ + str(hash(tuple(self.args)))

    def __repr__(self):
        return 'Key({!r}, *{!r})'.format(self.cls, self.args)

    def create(self):
        try:
            return self.cls(*self.args)
        except TypeError:
            raise ValueError('cannot create {!r} by {!r}'.format(self.cls, self.args))

    @property
    def descriptor_key(self):
        return self


class Descriptor(with_metaclass(ABCMeta, object)):
    '''
    abstruct base class of descriptors
    '''

    explicit_hydrogens = True
    gasteiger_charges = False
    kekulize = False

    @classmethod
    def preset(cls):
        '''
        generate preset descriptor instances

        Returns:
            iterable: descriptors
        '''
        yield cls()

    @classmethod
    def make_key(cls, *args):
        return Key(cls, *args)

    @property
    def dependencies(self):
        return None

    @property
    def descriptor_key(self):
        return self.make_key()

    @property
    def descriptor_name(self):
        '''
        get descriptor name

        Returns:
            str: descriptor name
        '''
        return str(self.descriptor_key)

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

    def _calculate(self, desc, cache):
        if desc.descriptor_key in cache:
            return cache[desc.descriptor_key]

        if isinstance(desc, Key):
            desc = desc.create()

        args = {name: self._calculate(dep, cache) if dep is not None else None
                for name, dep in (desc.dependencies or {}).items()}

        mol = self.molecule.get(
            explicitH=desc.explicit_hydrogens,
            gasteiger=desc.gasteiger_charges,
            kekulize=desc.kekulize,
        )
        r = desc.calculate(mol, **args)

        if desc.descriptor_key is None:
            raise ValueError('[bug] descriptor key not provided: {!r}'.format(desc))

        cache[desc.descriptor_key] = r
        return r

    def __call__(self, mol):
        '''
        calculate descriptors

        Returns:
            iterator<str, scalar>: iterator of descriptor name and value
        '''
        cache = {}
        self.molecule = Molecule(mol)

        return (
            (desc.descriptor_name, self._calculate(desc, cache))
            for desc in self.descriptors
        )


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
