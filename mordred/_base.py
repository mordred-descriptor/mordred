import math
import os
import sys

from abc import ABCMeta, abstractmethod
from importlib import import_module
from inspect import getsourcelines, isabstract
from sys import maxsize
from types import ModuleType

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
    r"""abstract base class of descriptors."""

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
        return '{}({})'.format(cls.__name__, ', '.join(map(pretty, args)))

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

    rtype = type(None)

    @classmethod
    def preset(cls):
        r"""generate preset descriptor instances.

        (abstract classmethod)

        :rtype: iterable
        """
        pass

    def dependencies(self):
        r"""descriptor dependencies.

        :rtype: {:py:class:`str`: (:py:class:`Descriptor` or :py:class:`None`)} or :py:class:`None`
        """
        pass

    @abstractmethod
    def calculate(self, mol):
        r"""calculate descriptor value.

        (abstract method)
        """
        pass

    def __call__(self, mol):
        r"""calculate single descriptor value.

        :returns: descriptor result
        :rtype: scalar
        """
        return Calculator(self)(mol)[0][1]

    @classmethod
    def is_descriptor_class(cls, desc):
        r"""check calculatable descriptor class or not.

        :rtype: :py:class:`bool`
        """
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

    :param descs: see :py:meth:`register` method
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

        :type descs: :py:class:`module`,
            :py:class:`Descriptor` class/instance or
            :py:class:`Iterable`

        :param descs: descriptors to register

            * :py:class:`module`: Descriptors in module
            * :py:class:`Descriptor` class: use :py:meth:`Descriptor.preset`

        """
        for desc in descs:
            if not hasattr(desc, '__iter__'):
                if Descriptor.is_descriptor_class(desc):
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
            cache[desc] = float('nan')
            return float('nan')

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

    def __call__(self, mol, error_callback=None):
        r"""calculate descriptors.

        :type mol: rdkit.Chem.Mol
        :param mol: molecular

        :type error_callback: callable
        :param error_callback: call when ransed Exception

        :rtype: [(Descriptor, scalar or nan)]
        :returns: iterator of descriptor and value
        """
        cache = {}
        self.molecule = Molecule(mol)

        if error_callback is None:
            def raise_error(e):
                raise e

            error_callback = raise_error

        rs = []
        for desc in self.descriptors:
            try:
                r = self._calculate(desc, cache)
            except Exception as e:
                r = error_callback(e)

            if not math.isnan(r) and not isinstance(r, desc.rtype):
                r = error_callback(DescriptorException(
                    desc,
                    TypeError('{!r}({}) is not {!r}'.format(r, type(r), desc.rtype)),
                    mol
                ))

            if math.isnan(r):
                rs.append((desc, float('nan')))
                continue

            rs.append((desc, desc.rtype(r)))

        return rs

    def _parallel(self, mols, processes, error_mode, callback, error_callback):
        from multiprocessing import Pool

        try:
            pool = Pool(
                processes,
                initializer=initializer,
                initargs=(self, error_mode),
            )

            kws = dict()

            if callback is not None:
                kws['callback'] = callback

            if error_callback is not None:
                kws['error_callback'] = error_callback

            def do_task(m):
                return pool.apply_async(
                    worker,
                    (m.ToBinary(),),
                    **kws
                )

            for m, result in [(m, do_task(m)) for m in mols]:

                if six.PY3:
                    yield m, result.get()
                else:
                    # timeout: avoid python2 KeyboardInterrupt bug.
                    # http://stackoverflow.com/a/1408476
                    yield m, result.get(1e9)

        finally:
            pool.terminate()
            pool.join()

    def _serial(self, mols, error_mode, callback, error_callback):
        calculate = make_calculator(self, error_mode)

        for m in mols:
            if error_callback is not None:
                try:
                    r = calculate(m)
                except Exception as e:
                    r = error_callback(e)
            else:
                r = calculate(m)

            if callback is not None:
                callback(r)

            yield m, r

    def map(self, mols, processes=None, error_mode='raise', callback=None, error_callback=None):
        r"""calculate descriptors over mols.

        :type mols: :py:class:`Iterable` (:py:class:`Mol`)
        :param mols: moleculars

        :type processes: :py:class:`int` or :py:class:`None`
        :param processes: number of process. None is :py:func:`multiprocessing.cpu_count`

        :type error_mode: :py:class:`str`
        :param error_mode:

            * 'raise': raise Exception
            * 'ignore': ignore Exception
            * 'log': print Exception to stderr and ingore Exception

        :type callback: :py:class:`Callable` ([(:py:class:`Descriptor`, scalar)])
            -> :py:class:`None`

        :param callback: call when calculate finished par molecule

        :type error_callback: :py:class:`Callable` (:py:class:`Exception`) -> scalar
        :param error_callback: call when Exception raised

        :rtype: :py:class:`Iterator` ((:py:class:`Mol`, [(:py:class:`Descriptor`, scalar)]]))
        """
        assert error_mode in set(['raise', 'ignore', 'log'])

        if processes == 1:
            return self._serial(mols, error_mode, callback, error_callback)
        else:
            return self._parallel(mols, processes, error_mode, callback, error_callback)


calculate = None


def initializer(calc, e_mode):
    global calculate

    calculate = make_calculator(calc, e_mode)


def make_calculator(calc, e_mode):
    if e_mode == 'raise':
        return calc

    elif e_mode == 'ignore':
        def ignore(e):
            return float('nan')

        return lambda m: calc(m, error_callback=ignore)

    else:
        def ignore_and_log(e):
            sys.stderr.write('{}\n'.format(e))
            return float('nan')

        return lambda m: calc(m, error_callback=ignore_and_log)


def worker(binary):
    return calculate(Chem.Mol(binary))


def all_descriptors():
    r"""yield all descriptors.

    :returns: all modules
    :rtype: :py:class:`Iterator` (:py:class:`module`)
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
    :param mdl: module to search

    :rtype: [:py:class:`Descriptor`]
    """
    descs = []

    for name in dir(mdl):
        if name[:1] == '_':
            continue

        desc = getattr(mdl, name)
        if Descriptor.is_descriptor_class(desc):
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
