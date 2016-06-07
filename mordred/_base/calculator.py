from .descriptor import Descriptor, NAException
from ..error import NA, Error, MultipleFragments
from types import ModuleType
from .context import Context
from inspect import getsourcelines
import sys
from .._util import Capture, DummyBar, NotebookWrapper
from itertools import chain
from contextlib import contextmanager
from tqdm import tqdm
from logging import getLogger


logger = getLogger('mordred')


class Calculator(object):
    r"""descriptor calculator.

    :param descs: see :py:meth:`register` method
    """

    __slots__ = (
        '_descriptors', '_explicit_hydrogens', '_kekulizes', '_require_3D',
        '_cache', '_progress_bar'
    )

    def __setstate__(self, dict):
        self._descriptors = dict.get('_descriptors', [])
        self._explicit_hydrogens = dict.get('_explicit_hydrogens', set([True, False]))
        self._kekulizes = dict.get('_kekulizes', set([True, False]))
        self._require_3D = dict.get('_require_3D', False)

    def __reduce_ex__(self, version):
        return self.__class__, (), {
            '_descriptors': self._descriptors,
            '_explicit_hydrogens': self._explicit_hydrogens,
            '_kekulizes': self._kekulizes,
            '_require_3D': self._require_3D,
        }

    def __init__(self, descs=[], exclude3D=False):
        self._descriptors = []

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
        self._descriptors = []
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

    def _calculate_one(self, cxt, desc, reset):
        if desc in self._cache:
            return self._cache[desc]

        if reset:
            cxt.reset()
        desc._context = cxt
        cxt.add_stack(desc)

        if desc.require_connected and desc._context.n_frags != 1:
            desc.fail(MultipleFragments())

        args = {
            name: self._calculate_one(cxt, dep, False)
            if dep is not None else None
            for name, dep in (desc.dependencies() or {}).items()
        }

        r = desc.calculate(**args)

        self._check_rtype(desc, r)

        self._cache[desc] = r

        return r

    def _check_rtype(self, desc, result):
        if desc.rtype is None:
            return

        if isinstance(result, Error):
            return

        if not isinstance(result, desc.rtype):
            pass  # TODO

    def _calculate(self, cxt):
        self._cache = {}
        for desc in self.descriptors:
            try:
                yield self._calculate_one(cxt, desc, True)
            except NAException as e:
                yield NA(e.error, desc._context.get_stack())
            except Exception as e:
                yield Error(e, desc._context.get_stack())

    def __call__(self, mol, id=-1):
        r"""calculate descriptors.

        :type mol: rdkit.Chem.Mol
        :param mol: molecular

        :type id: int
        :param id: conformer id

        :rtype: [scalar or Error]
        :returns: iterator of descriptor and value
        """
        return list(self._calculate(Context.from_calculator(self, mol, id)))

    def _serial(self, mols, nmols, quiet, ipynb, id):
        with self._progress(quiet, nmols, ipynb) as bar:
            for m in mols:
                with Capture() as capture:
                    r = list(self._calculate(Context.from_calculator(self, m, id)))

                for e in capture.result:
                    e = e.rstrip()
                    if not e:
                        continue

                    bar.write(e, file=capture.orig)

                yield m, r
                bar.update()

    @contextmanager
    def _progress(self, quiet, total, ipynb):
        args = {
            'dynamic_ncols': True,
            'leave': True,
            'total': total
        }

        if quiet:
            Bar = DummyBar
        elif ipynb:
            Bar = NotebookWrapper
        else:
            Bar = tqdm

        try:
            with Bar(**args) as self._progress_bar:
                yield self._progress_bar
        finally:
            del self._progress_bar

    def logging(self, s):
        p = getattr(self, '_progress_bar', None)
        if p is not None:
            p.write(s, file=sys.stderr)
            return

        logger.warn(s)

    def map(self, mols, nproc=None, nmols=None, quiet=False, ipynb=False, id=-1):
        r"""calculate descriptors over mols.

        :type mols: :py:class:`Iterable` (:py:class:`Mol`)
        :param mols: moleculars

        :type nproc: :py:class:`int` or :py:class:`None`
        :param nproc: number of process. None is :py:func:`multiprocessing.cpu_count`

        :type nmols: :py:class:`None` or :py:class:`int`
        :param nmols: number of all mols for display progress bar

        :type quiet: :py:class:`bool`
        :param quiet: suppress progress bar

        :type ipynb: :py:class:`bool`
        :param ipynb: use ipython notebook progress bar

        :type id: :py:class:`int`
        :param id: conformer id

        :rtype: :py:class:`Iterator` ((:py:class:`Mol`, [scalar]]))
        """

        if hasattr(mols, '__len__'):
            nmols = len(mols)

        if nproc == 1:
            return self._serial(mols, nmols=nmols, quiet=quiet, ipynb=ipynb, id=id)
        else:
            return self._parallel(mols, nproc, nmols=nmols, quiet=quiet, ipynb=ipynb, id=id)

    def pandas(self, mols, mol_name='mol',
               nproc=None, nmols=None, quiet=False, ipynb=False, id=-1):
        r"""calculate descriptors over mols.

        :type mol_name: str
        :param mol_name: molecular column name

        :rtype: :py:class:`pandas.DataFrame`
        """
        import pandas

        return pandas.DataFrame(
            (chain([m], d) for m, d in self.map(mols, nproc, nmols, quiet, ipynb, id)),
            columns=list(chain([mol_name], (str(d) for d in self.descriptors)))
        )


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
            return sys.maxsize

    descs.sort(key=key_by_def)
    return descs
