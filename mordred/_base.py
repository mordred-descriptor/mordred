import io
import os
import sys

from abc import ABCMeta, abstractmethod
from importlib import import_module
from inspect import getsourcelines, isabstract
from sys import maxsize
from types import ModuleType

import numpy as np

from rdkit import Chem
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges

import six


class MordredException(Exception):
    pass


class MordredAttributeError(AttributeError, MordredException):
    def __init__(self, desc, args):
        super(AttributeError, self).__init__()
        self.desc = desc
        self.args = args

    def __reduce_ex__(self, version):
        return self.__class__, (self.desc, self.args)

    def __str__(self):
        return '{}({})'.format(self.args, self.desc)


class DescriptorException(MordredException):
    def __init__(self, desc, e, mol, parent=None):
        self.desc = desc
        self.e = e
        self.mol = mol
        self.parent = parent

    def __reduce_ex__(self, version):
        return self.__class__, (self.desc, self.e, self.mol, self.parent)

    def __str__(self):
        if self.parent is None:
            return '{}({!r}): {}'.format(
                self.desc,
                Chem.MolToSmiles(self.mol),
                self.e,
            )

        return '{}/{}({!r}): {}'.format(
            self.parent,
            self.desc,
            Chem.MolToSmiles(self.mol),
            self.e,
        )


def pretty(a):
    p = getattr(a, 'name', None)
    return repr(a if p is None else p)


class Descriptor(six.with_metaclass(ABCMeta, object)):
    r"""abstruct base class of descriptors."""

    explicit_hydrogens = True
    gasteiger_charges = False
    kekulize = False
    require_connected = False

    _reduce_ex_version = 3

    @abstractmethod
    def __reduce_ex__(self, version):
        pass

    def __repr__(self):
        cls, args = self.__reduce_ex__(self._reduce_ex_version)
        return '{}({})'.format(cls, ', '.join(map(pretty, args)))

    def __hash__(self):
        return hash(self.__reduce_ex__(self._reduce_ex_version))

    def __eq__(self, other):
        l = self.__reduce_ex__(self._reduce_ex_version)
        r = other.__reduce_ex__(self._reduce_ex_version)
        return l.__eq__(r)

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        l = self.__reduce_ex__(self._reduce_ex_version)
        r = other.__reduce_ex__(self._reduce_ex_version)
        return l.__lt__(r)

    @classmethod
    def preset(cls):
        r"""generate preset descriptor instances.

        (abstruct classmethod)

        :rtype: iterable
        """
        pass

    def dependencies(self):
        r"""descriptor dependencies.

        :rtype: {str: (Descriptor or None)} or None
        """
        pass

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
        return Calculator(self)(mol)[0][1]

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
            raise DescriptorException(desc, e, mol, parent)

        cache[desc] = r
        return r

    def __call__(self, mol, on_exception='raise'):
        r"""calculate descriptors.

        :type mol: rdkit.Chem.Mol
        :param mol: molecular

        :type on_exception: :py:class:`str`, :py:class:`io.TextIOWrapper` or callable
        :param on_exception:

            * 'raise': raise Exception
            * 'ignore': ignore error
            * 'log_and_ignore': log exception to stderr and ignore error
            * io.TextIOWrapper: log exception to file and ignore error
            * callable: call with exception

        :rtype: [(Descriptor, scalar or nan)]
        :returns: iterator of descriptor and value
        """
        cache = {}
        self.molecule = Molecule(mol)

        if on_exception == 'raise':
            def handler(e):
                raise e
        elif on_exception == 'ignore':
            def handler(e):
                pass
        elif on_exception == 'log_and_ignore':
            def handler(e):
                sys.stderr.write('{}\n'.format(e))
        elif isinstance(on_exception, io.TextIOWrapper) and on_exception.writable():
            def handler(e):
                on_exception.write('{}\n'.format(e))
        else:
            def handler(e):
                on_exception(e)

        rs = []
        for desc in self.descriptors:
            try:
                r = self._calculate(desc, cache)
            except Exception as e:
                handler(e)
                r = np.nan

            if not isinstance(
                    r,
                    (six.integer_types, np.integer,
                     float, np.floating,
                     bool, np.bool_)):

                handler(DescriptorException(
                    desc,
                    ValueError('not int or float: {!r}({})'.format(r, type(r))),
                    mol
                ))

                r = np.nan

            rs.append((desc, r))

        return rs

    def _parallel(self, mols, processes=None, on_exception='raise'):
        from multiprocessing import Pool

        try:
            pool = Pool(
                processes,
                initializer=initializer,
                initargs=(self, on_exception),
            )

            for m, result in [
                    (m, pool.apply_async(worker, (m.ToBinary(),)))
                    for m in mols]:

                if six.PY3:
                    yield m, result.get()
                else:
                    # timeout: avoid python2 KeyboardInterrupt bug.
                    # http://stackoverflow.com/a/1408476
                    yield m, result.get(1e9)

        finally:
            pool.terminate()
            pool.join()

    def map(self, mols, processes=None, on_exception='raise'):
        r"""calculate descriptors over mols.

        :param mols: moleculars
        :type mols: iterable(rdkit.Chem.Mol)

        :param processes: number of process. None is multiprocessing.cpu_count()
        :type processes: int or None

        :rtype: iterator((rdkit.Chem.Mol, [(Descriptor, scalar)]]))
        """
        if processes == 1:
            return ((m, list(self(m, on_exception=on_exception))) for m in mols)
        else:
            return self._parallel(mols, processes, on_exception)


calculate = None


def initializer(calc, on_exception):
    global calculate
    calculate = lambda m: calc(m, on_exception=on_exception)


def worker(binary):
    return calculate(Chem.Mol(binary))


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
