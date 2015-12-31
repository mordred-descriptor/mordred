from rdkit import Chem
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges
from six import with_metaclass
from abc import ABCMeta, abstractproperty, abstractmethod
from types import ModuleType

import inspect


__all__ =\
    'Calculator',\
    'Descriptor',


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
        if len(self.args) == 0:
            return self.cls.__name__
        else:
            return self.cls.__name__ + str(hash(self.args))

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


class DescriptorMeta(ABCMeta):
    def __new__(cls, name, base, dic):
        cls = ABCMeta.__new__(cls, name, base, dic)
        if 'descriptor_defaults' not in dic:
            spec = inspect.getfullargspec(cls.__init__)
            if len(spec.args) == 1:
                cls.descriptor_defaults = [()]
            elif len(spec.args) - 1 == len(spec.defaults or []):
                cls.descriptor_defaults = [spec.defaults]

        return cls


class Descriptor(with_metaclass(DescriptorMeta, object)):
    explicitHydrogens = True
    gasteigerCharges = False
    kekulize = False

    descriptor_defaults = []

    @classmethod
    def make_key(cls, *args):
        return Key(cls, *args)

    @property
    def dependencies(self):
        return dict()

    @abstractproperty
    def descriptor_key(self):
        pass

    @property
    def descriptor_name(self):
        return str(self.descriptor_key)

    @abstractmethod
    def calculate(self, mol):
        pass

    def __call__(self, mol):
        return Calculator(self)(mol)[0][1]


class Molecule(object):
    def __init__(self, orig):
        Chem.Kekulize(orig)
        self.orig = orig
        self.cache = dict()

    def key(self, explicitH=True, kekulize=False, gasteiger=False):
        return (explicitH, kekulize, gasteiger)

    def hydration(self, explicitH):
        key = self.key(explicitH=explicitH)
        if key in self.cache:
            return self.cache[key]

        mol = Chem.AddHs(self.orig) if explicitH else Chem.RemoveHs(self.orig)
        self.cache[key] = mol
        return mol

    def kekulize(self, explicitH):
        key = self.key(explicitH=explicitH, kekulize=True)
        if key in self.cache:
            return self.cache[key]

        mol = Chem.Mol(self.hydration(explicitH))
        Chem.Kekulize(mol)
        self.cache[key] = mol
        return mol

    def gasteiger(self, explicitH, kekulize):
        key = self.key(explicitH=explicitH, kekulize=kekulize, gasteiger=True)
        if key in self.cache:
            return self.cache[key]

        mol = self.kekulize(explicitH) if kekulize else self.hydration(explicitH)
        ComputeGasteigerCharges(mol)
        self.cache[key] = mol
        return mol

    def get(self, explicitH, kekulize, gasteiger):
        key = self.key(explicitH=explicitH, kekulize=kekulize, gasteiger=gasteiger)
        if key in self.cache:
            return self.cache[key]

        if gasteiger:
            return self.gasteiger(explicitH, kekulize)
        elif kekulize:
            return self.kekulize(explicitH)
        else:
            return self.hydration(explicitH)


class Calculator(object):
    def __init__(self, *descs):
        self.descriptors = []
        self.explicitH = False
        self.gasteiger = False
        self.kekulize = False

        self.register(*descs)

    def _register_one(self, desc):
        # check desc is descriptor or not
        if not isinstance(desc, Descriptor):
            raise ValueError('{!r} is not descriptor'.format(desc))

        self.descriptors.append(desc)

        if desc.explicitHydrogens:
            self.explicitH = True

        if desc.gasteigerCharges:
            self.gasteiger = True

        if desc.kekulize:
            self.kekulize = True

    def register(self, *descs):
        for desc in descs:
            if not hasattr(desc, '__iter__'):
                if isinstance(desc, type):
                    for args in desc.descriptor_defaults:
                        self._register_one(desc(*args))

                elif isinstance(desc, ModuleType):
                    for name in dir(desc):
                        if name[:1] == '_':
                            continue

                        d = getattr(desc, name)
                        if issubclass(d, Descriptor):
                            self.register(d)

                else:
                    self._register_one(desc)

            elif isinstance(desc, tuple) and isinstance(desc[0], type):
                self._register_one(desc[0](*desc[1:]))

            else:
                for d in desc:
                    self.register(d)

    def _calculate(self, desc, cache):
        if desc.descriptor_key in cache:
            return cache[desc.descriptor_key]

        if isinstance(desc, Key):
            desc = desc.create()

        args = {name: self._calculate(dep, cache)
                for name, dep in desc.dependencies.items()}

        mol = self.molecule.get(
            explicitH=desc.explicitHydrogens,
            gasteiger=desc.gasteigerCharges,
            kekulize=desc.kekulize,
        )
        r = desc.calculate(mol, **args)

        if desc.descriptor_key is None:
            raise ValueError('[bug] descriptor key not provided: {!r}'.format(desc))

        cache[desc.descriptor_key] = r
        return r

    def __call__(self, mol):
        cache = {}
        self.molecule = Molecule(mol)

        for desc in self.descriptors:
            self._calculate(desc, cache)

        return [(d.descriptor_name, cache[d.descriptor_key]) for d in self.descriptors]
