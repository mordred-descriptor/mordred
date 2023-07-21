from __future__ import print_function

import sys
import warnings
from types import ModuleType
from contextlib import contextmanager
from multiprocessing import cpu_count
from packaging.version import Version as StrictVersion

from .._util import Capture, DummyBar
from ..error import Error, Missing, MultipleFragments, DuplicatedDescriptorName
from .result import Result
from .context import Context
from .descriptor import Descriptor, MissingValueException, is_descriptor_class

try:
    from importlib.metadata import version as importlib_version
except ImportError:
    from importlib_metadata import version as importlib_version

__version__ = importlib_version("mordredcommunity")

try:
    from tqdm import tqdm
    from .._util import NotebookWrapper
except ImportError:
    tqdm = NotebookWrapper = DummyBar


class Calculator(object):
    r"""descriptor calculator.

    Parameters:
        descs: see Calculator.register() method
        ignore_3D: see Calculator.register() method

    """

    __slots__ = (
        "_descriptors",
        "_name_dict",
        "_explicit_hydrogens",
        "_kekulizes",
        "_require_3D",
        "_cache",
        "_debug",
        "_progress_bar",
        "_config",
    )

    def __setstate__(self, dict):
        ds = self._descriptors = dict.get("_descriptors", [])
        self._name_dict = {str(d): d for d in ds}
        self._explicit_hydrogens = dict.get("_explicit_hydrogens", {True, False})
        self._kekulizes = dict.get("_kekulizes", {True, False})
        self._require_3D = dict.get("_require_3D", False)

    @classmethod
    def from_json(cls, obj):
        """Create Calculator from json descriptor objects.

        Parameters:
            obj(list or dict): descriptors to register

        Returns:
            Calculator: calculator

        """
        calc = cls()
        calc.register_json(obj)
        return calc

    def register_json(self, obj):
        """Register Descriptors from json descriptor objects.

        Parameters:
            obj(list or dict): descriptors to register

        """
        if not isinstance(obj, list):
            obj = [obj]

        self.register(Descriptor.from_json(j) for j in obj)

    def to_json(self):
        """Convert descriptors to json serializable data.

        Returns:
            list: descriptors

        """
        return [d.to_json() for d in self.descriptors]

    def __reduce_ex__(self, version):
        return (
            self.__class__,
            (),
            {
                "_config": self._config,
                "_descriptors": self._descriptors,
                "_explicit_hydrogens": self._explicit_hydrogens,
                "_kekulizes": self._kekulizes,
                "_require_3D": self._require_3D,
            },
        )

    def __getitem__(self, key):
        return self._name_dict[key]

    def __init__(self, descs=None, version=None, ignore_3D=False, config=None):
        if descs is None:
            descs = []

        if config is None:
            config = {}

        self._descriptors = []
        self._name_dict = {}

        self._explicit_hydrogens = set()
        self._kekulizes = set()
        self._require_3D = False
        self._debug = False
        self._config = config

        self.register(descs, version=version, ignore_3D=ignore_3D)

    def config(self, **configs):
        r"""Set global configuration."""
        self._config.update(configs)

    @property
    def descriptors(self):
        r"""All descriptors.

        you can get/set/delete descriptor.

        Returns:
            tuple[Descriptor]: registered descriptors

        """
        return tuple(self._descriptors)

    @descriptors.setter
    def descriptors(self, descs):
        del self.descriptors
        self.register(descs)

    @descriptors.deleter
    def descriptors(self):
        self._descriptors = []
        self._name_dict = {}
        self._explicit_hydrogens.clear()
        self._kekulizes.clear()
        self._require_3D = False

    def __len__(self):
        return len(self._descriptors)

    def _register_one(self, desc, check_only=False, ignore_3D=False):
        if not isinstance(desc, Descriptor):
            raise ValueError("{!r} is not descriptor".format(desc))

        if ignore_3D and desc.require_3D:
            return

        self._explicit_hydrogens.add(bool(desc.explicit_hydrogens))
        self._kekulizes.add(bool(desc.kekulize))
        self._require_3D |= desc.require_3D

        for dep in (desc.dependencies() or {}).values():
            if isinstance(dep, Descriptor):
                self._register_one(dep, check_only=True)

        if not check_only:
            sdesc = str(desc)
            old = self._name_dict.get(sdesc)
            if old is not None:
                raise DuplicatedDescriptorName(desc, old)

            self._name_dict[sdesc] = desc
            self._descriptors.append(desc)

    def register(self, desc, version=None, ignore_3D=False):
        r"""Register descriptors.

        Descriptor-like:
            * Descriptor instance: self
            * Descriptor class: use Descriptor.preset() method
            * module: use Descriptor-likes in module
            * Iterable: use Descriptor-likes in Iterable

        Parameters:
            desc(Descriptor-like): descriptors to register
            version(str): version
            ignore_3D(bool): ignore 3D descriptors

        """
        if version is None:
            version = __version__

        version = StrictVersion(version)

        return self._register(desc, version, ignore_3D)

    def _register(self, desc, version, ignore_3D):
        if not hasattr(desc, "__iter__"):
            if is_descriptor_class(desc):
                if desc.since > version:
                    return

                for d in desc.preset(version=version):
                    self._register_one(d, ignore_3D=ignore_3D)

            elif isinstance(desc, ModuleType):
                self._register(
                    get_descriptors_in_module(desc),
                    version=version,
                    ignore_3D=ignore_3D,
                )

            else:
                self._register_one(desc, ignore_3D=ignore_3D)

        else:
            for d in desc:
                self._register(d, version=version, ignore_3D=ignore_3D)

    def _calculate_one(self, cxt, desc, reset):
        if desc in self._cache:
            return self._cache[desc]

        if reset:
            cxt.reset()
        desc._context = cxt
        cxt.add_stack(desc)

        if desc.require_connected and desc._context.n_frags != 1:
            return False, Missing(MultipleFragments(), desc._context.get_stack())

        args = {}
        for name, dep in (desc.dependencies() or {}).items():
            if dep is None:
                args[name] = None
            else:
                ok, r = self._calculate_one(cxt, dep, False)
                if ok:
                    args[name] = r
                else:
                    return False, r

        ok = False
        try:
            r = desc.calculate(**args)
            if self._debug:
                self._check_rtype(desc, r)
            ok = True
        except MissingValueException as e:
            r = Missing(e.error, desc._context.get_stack())
        except Exception as e:
            r = Error(e, desc._context.get_stack())

        self._cache[desc] = ok, r

        return ok, r

    def _check_rtype(self, desc, result):
        if desc.rtype is None:
            return

        if isinstance(result, Error):
            return

        if not isinstance(result, desc.rtype):
            raise TypeError("{} not match {}".format(result, desc.rtype))

    def _calculate(self, cxt):
        self._cache = {}
        for desc in self.descriptors:
            _, r = self._calculate_one(cxt, desc, True)
            yield r

    def __call__(self, mol, id=-1):
        r"""Calculate descriptors.

        :type mol: rdkit.Chem.Mol
        :param mol: molecular

        :type id: int
        :param id: conformer id

        :rtype: Result[scalar or Error]
        :returns: iterator of descriptor and value
        """
        return self._wrap_result(
            mol, self._calculate(Context.from_calculator(self, mol, id))
        )

    def _wrap_result(self, mol, r):
        return Result(mol, r, self._descriptors)

    def _serial(self, mols, nmols, quiet, ipynb, id):
        with self._progress(quiet, nmols, ipynb) as bar:
            for m in mols:
                with Capture() as capture:
                    r = self._wrap_result(
                        m, self._calculate(Context.from_calculator(self, m, id))
                    )

                for e in capture.result:
                    e = e.rstrip()
                    if not e:
                        continue

                    bar.write(e, file=capture.orig)

                yield r
                bar.update()

    @contextmanager
    def _progress(self, quiet, total, ipynb):
        args = {"dynamic_ncols": True, "leave": True, "total": total}

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
            if hasattr(self, "_progress_bar"):
                del self._progress_bar

    def echo(self, s, file=sys.stdout, end="\n"):
        """Output message.

        Parameters:
            s(str): message to output
            file(file-like): output to
            end(str): end mark of message

        Return:
            None

        """
        p = getattr(self, "_progress_bar", None)
        if p is not None:
            p.write(s, file=file, end="\n")
            return

        print(s, file=file, end="\n")  # noqa: T003

    def map(self, mols, nproc=None, nmols=None, quiet=False, ipynb=False, id=-1):
        r"""Calculate descriptors over mols.

        Parameters:
            mols(Iterable[rdkit.Mol]): moleculars

            nproc(int): number of process to use. default: multiprocessing.cpu_count()

            nmols(int): number of all mols to use in progress-bar. default: mols.__len__()

            quiet(bool): don't show progress bar. default: False

            ipynb(bool): use ipython notebook progress bar. default: False

            id(int): conformer id to use. default: -1.

        Returns:
            Iterator[Result[scalar]]

        """
        if nproc is None:
            nproc = cpu_count()

        if hasattr(mols, "__len__"):
            nmols = len(mols)

        if nproc == 1:
            return self._serial(mols, nmols=nmols, quiet=quiet, ipynb=ipynb, id=id)
        else:
            return self._parallel(
                mols, nproc, nmols=nmols, quiet=quiet, ipynb=ipynb, id=id
            )

    def pandas(self, mols, nproc=None, nmols=None, quiet=False, ipynb=False, id=-1):
        r"""Calculate descriptors over mols.

        Returns:
            pandas.DataFrame

        """
        from .pandas_module import MordredDataFrame, Series

        if isinstance(mols, Series):
            index = mols.index
        else:
            index = None

        return MordredDataFrame(
            (list(r) for r in self.map(mols, nproc, nmols, quiet, ipynb, id)),
            columns=[str(d) for d in self.descriptors],
            index=index,
        )


def get_descriptors_in_module(mdl, submodule=True):
    r"""Get descriptors in module.

    Parameters:
        mdl(module): module to search
        submodule(bool): search recursively

    Returns:
        Iterator[Descriptor]

    """
    __all__ = getattr(mdl, "__all__", None)
    if __all__ is None:
        __all__ = dir(mdl)

    all_values = (getattr(mdl, name) for name in __all__ if name[:1] != "_")

    if submodule:
        for v in all_values:
            if is_descriptor_class(v):
                yield v
            if isinstance(v, ModuleType):
                for v in get_descriptors_in_module(v, submodule=True):
                    yield v

    else:
        for v in all_values:
            if is_descriptor_class(v):
                yield v
