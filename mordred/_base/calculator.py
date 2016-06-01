from .descriptor import Descriptor
from logging import getLogger
from types import ModuleType
from .exception import FragmentError
from numbers import Real
import sys
import os
import traceback
from .context import Context
from inspect import getsourcelines
from sys import maxsize


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

    def __init__(self, descs=[], exclude3D=False):
        self._descriptors = []
        self.logger = getLogger(__name__)

        self._explicit_hydrogens = set()
        self._kekulizes = set()
        self._require_3D = False

        self.register(descs, exclude3D=exclude3D)

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

    def _register_one(self, desc, check_only=False, exclude3D=False):
        if not isinstance(desc, Descriptor):
            raise ValueError('{!r} is not descriptor'.format(desc))

        if exclude3D and desc.require_3D:
            return

        self._explicit_hydrogens.add(bool(desc.explicit_hydrogens))
        self._kekulizes.add(bool(desc.kekulize))
        self._require_3D |= desc.require_3D

        for dep in (desc.dependencies() or {}).values():
            if isinstance(dep, Descriptor):
                self._register_one(dep, check_only=True)

        if not check_only:
            self._descriptors.append(desc)

    def register(self, desc, exclude3D=False):
        r"""register descriptors.

        :type desc: :py:class:`module`,
            :py:class:`Descriptor` class/instance or
            :py:class:`Iterable`

        :param desc: descriptors to register

            * :py:class:`module`: Descriptors in module
            * :py:class:`Descriptor` class: use :py:meth:`Descriptor.preset`
        """
        if not hasattr(desc, '__iter__'):
            if Descriptor.is_descriptor_class(desc):
                for d in desc.preset():
                    self._register_one(d, exclude3D=exclude3D)

            elif isinstance(desc, ModuleType):
                self.register(get_descriptors_from_module(desc), exclude3D=exclude3D)

            else:
                self._register_one(desc, exclude3D=exclude3D)

        else:
            for d in desc:
                self.register(d, exclude3D=exclude3D)

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
