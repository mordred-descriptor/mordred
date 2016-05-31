import os
import sys

from abc import ABCMeta, abstractmethod
from importlib import import_module
from inspect import getsourcelines, isabstract
from sys import maxsize
from types import ModuleType
from numbers import Real

from rdkit import Chem

import traceback
import six
from logging import getLogger

from ._exception import FragmentError
from ._util import conformer_to_numpy


class Descriptor(six.with_metaclass(ABCMeta, object)):
    r"""abstract base class of descriptors."""

    explicit_hydrogens = True
    kekulize = False
    require_connected = False
    require_3D = False

    _reduce_ex_version = 3

    @abstractmethod
    def __reduce_ex__(self, version):
        pass

    @staticmethod
    def _pretty(a):
        p = getattr(a, 'name', None)
        return repr(a if p is None else p)

    def __repr__(self):
        cls, args = self.__reduce_ex__(self._reduce_ex_version)
        return '{}({})'.format(cls.__name__, ', '.join(map(self._pretty, args)))

    def __hash__(self):
        return hash(self.__reduce_ex__(self._reduce_ex_version))

    def __compare_by_reduce(meth):
        def compare(self, other):
            l = self.__reduce_ex__(self._reduce_ex_version)
            r = other.__reduce_ex__(other._reduce_ex_version)
            return getattr(l, meth)(r)

        return compare

    __eq__ = __compare_by_reduce('__eq__')
    __ne__ = __compare_by_reduce('__ne__')

    __lt__ = __compare_by_reduce('__lt__')
    __gt__ = __compare_by_reduce('__gt__')
    __le__ = __compare_by_reduce('__le__')
    __ge__ = __compare_by_reduce('__ge__')

    rtype = None

    @classmethod
    def preset(cls):
        r"""generate preset descriptor instances.

        :rtype: iterable
        """
        return ()

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

    def __call__(self, mol, coord_id=-1):
        r"""calculate single descriptor value.

        :returns: descriptor result
        :rtype: scalar
        """
        return Calculator(self)(mol, coord_id)[0]

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


class Context(object):
    def __reduce_ex__(self, version):
        return self.__class__, (None,), {
            '_mols': self._mols,
            '_coords': self._coords,
            'n_frags': self.n_frags,
            'name': self.name,
        }

    def __str__(self):
        return self.name

    @classmethod
    def from_calculator(cls, calc, mol, id):
        return cls(mol, calc._require_3D, calc._explicit_hydrogens, calc._kekulizes, id)

    __tf = set([True, False])

    def __init__(self, mol, require_3D=False, explicit_hydrogens=__tf, kekulizes=__tf, id=-1):
        if mol is None:
            return

        self._mols = {}
        self._coords = {}

        self.n_frags = len(Chem.GetMolFrags(mol))

        if mol.HasProp('_Name'):
            self.name = mol.GetProp('_Name')
        else:
            self.name = Chem.MolToSmiles(Chem.RemoveHs(mol))

        for eh, ke in ((eh, ke) for eh in explicit_hydrogens for ke in kekulizes):
            m = (Chem.AddHs if eh else Chem.RemoveHs)(mol)

            if ke:
                Chem.Kekulize(m)

            if require_3D:
                self._coords[eh, ke] = conformer_to_numpy(m.GetConformer(id))

            m.RemoveAllConformers()
            self._mols[eh, ke] = m

    def get_coord(self, desc):
        return self._coords[desc.explicit_hydrogens, desc.kekulize]

    def get_mol(self, desc):
        return self._mols[desc.explicit_hydrogens, desc.kekulize]


class Calculator(object):
    r"""descriptor calculator.

    :param descs: see :py:meth:`register` method
    """

    def __reduce_ex__(self, version):
        return self.__class__, (), {
            '_descriptors': self._descriptors,
            '_explicit_hydrogens': self._explicit_hydrogens,
            '_kekulizes': self._kekulizes,
            '_require_3D': self._require_3D,
        }

    def __init__(self, *descs):
        self._descriptors = []
        self.logger = getLogger(__name__)

        self._explicit_hydrogens = set()
        self._kekulizes = set()
        self._require_3D = False

        self.register(*descs)

    @property
    def descriptors(self):
        r'''all descriptors.

        you can get/set/delete descriptor.
        '''
        return tuple(self._descriptors)

    @descriptors.setter
    def descriptors(self, descs):
        del self.descriptors
        self.register(descs)

    @descriptors.deleter
    def descriptors(self):
        self._descriptors[:] = []
        self._explicit_hydrogens.clear()
        self._kekulizes.clear()
        self._require_3D = False

    def __len__(self):
        return len(self._descriptors)

    def _register_one(self, desc, check_only=False):
        if not isinstance(desc, Descriptor):
            raise ValueError('{!r} is not descriptor'.format(desc))

        self._explicit_hydrogens.add(bool(desc.explicit_hydrogens))
        self._kekulizes.add(bool(desc.kekulize))
        self._require_3D |= desc.require_3D

        for dep in (desc.dependencies() or {}).values():
            if isinstance(dep, Descriptor):
                self._register_one(dep, True)

        if not check_only:
            self._descriptors.append(desc)

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

    def _calculate_one(self, cxt, desc, caller=None):
        if desc in self._cache:
            return self._cache[desc]

        if caller is None:
            caller = desc

        if desc.require_connected and cxt.n_frags != 1:
            raise FragmentError(cxt, caller)

        args = {
            name: self._calculate_one(cxt, dep, caller)
            if dep is not None else None
            for name, dep in (desc.dependencies() or {}).items()
        }

        mol = cxt.get_mol(desc)

        if desc.require_3D:
            r = desc.calculate(mol, cxt.get_coord(desc), **args)
        else:
            r = desc.calculate(mol, **args)

        if desc.rtype is not None and (not isinstance(r, Real) or not isinstance(r, desc.rtype)):
            self.logger.debug(
                '%s excepted returning %s, but returns %s',
                repr(desc), desc.rtype.__name__, repr(r)
            )

        self._cache[desc] = r

        return r

    def _calculate(self, cxt):
        self._cache = {}
        self._exceptions = set()
        try:
            for desc in self.descriptors:
                try:
                    yield self._calculate_one(cxt, desc)
                except Exception as e:
                    if e.__class__ not in self._exceptions:
                        exc_type, exc_value, exc_traceback = sys.exc_info()
                        tbs = traceback.extract_tb(exc_traceback)[-1:]
                        for tb in tbs:
                            filename, line, _, text = tb
                            self.logger.warning(
                                '%s:%d %s',
                                os.path.basename(filename), line, str(e)
                            )
                        self._exceptions.add(e.__class__)

                    yield float('nan')

        finally:
            del self._cache
            del self._exceptions

    def __call__(self, mol, id=-1):
        r"""calculate descriptors.

        :type mol: rdkit.Chem.Mol
        :param mol: molecular

        :type id: int
        :param id: conformer id

        :rtype: [scalar or nan]
        :returns: iterator of descriptor and value
        """
        return list(self._calculate(Context.from_calculator(self, mol, id)))

    def _parallel(self, mols, processes, callback):
        from multiprocessing import Pool

        try:
            pool = Pool(
                processes,
                initializer=initializer,
                initargs=(self,),
            )

            if callback is None:
                def do_task(mol):
                    args = Context.from_calculator(self, mol, -1),
                    return pool.apply_async(worker, args)
            else:
                def do_task(mol):
                    args = Context.from_calculator(self, mol, -1),
                    return pool.apply_async(worker, args, callback=callback)

            if six.PY3:
                def get_result(r):
                    return r.get()
            else:
                # timeout: avoid python2 KeyboardInterrupt bug.
                # http://stackoverflow.com/a/1408476
                def get_result(r):
                    return r.get(1e9)

            for m, result in [(m, do_task(m)) for m in mols]:
                yield m, get_result(result)

        finally:
            pool.terminate()
            pool.join()

    def _serial(self, mols, callback):
        for m in mols:
            r = list(self._calculate(Context.from_calculator(self, m, -1)))
            if callback is not None:
                callback((m, r))

            yield m, r

    def map(self, mols, processes=None, callback=None):
        r"""calculate descriptors over mols.

        :type mols: :py:class:`Iterable` (:py:class:`Mol`)
        :param mols: moleculars

        :type processes: :py:class:`int` or :py:class:`None`
        :param processes: number of process. None is :py:func:`multiprocessing.cpu_count`

        :type callback: :py:class:`Callable` ([scalar])
            -> :py:class:`None`

        :param callback: call when calculate finished par molecule

        :rtype: :py:class:`Iterator` ((:py:class:`Mol`, [scalar]]))
        """
        if processes == 1:
            return self._serial(mols, callback=callback)
        else:
            return self._parallel(mols, processes, callback=callback)


calculator = None


def initializer(calc):
    global calculator

    calculator = calc


def worker(cxt):
    return list(calculator._calculate(cxt))


def all_modules():
    base_dir = os.path.dirname(__file__)

    for name in os.listdir(base_dir):
        name, ext = os.path.splitext(name)
        if name[:1] == '_' or ext != '.py':
            continue

        yield import_module('.' + name, __package__)


def all_descriptors(with_3D=True):
    r"""yield all descriptors.

    :returns: all modules
    :rtype: :py:class:`Iterator` (:py:class:`Descriptor`)
    """

    for mdl in all_modules():
        for desc in get_descriptors_from_module(mdl):
            if not with_3D and desc.require_3D is True:
                continue

            yield desc


def get_descriptors_from_module(mdl):
    r"""get descriptors from module.

    :type mdl: module
    :param mdl: module to search

    :rtype: [:py:class:`Descriptor`]
    """

    __all__ = getattr(mdl, '__all__', None)
    if __all__ is None:
        __all__ = dir(mdl)

    descs = [
        fn
        for fn in (getattr(mdl, name) for name in __all__ if name[:1] != '_')
        if Descriptor.is_descriptor_class(fn)
    ]

    def key_by_def(d):
        try:
            return getsourcelines(d)[1]
        except IOError:
            return maxsize

    descs.sort(key=key_by_def)
    return descs
