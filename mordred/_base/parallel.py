from itertools import islice
from collections import deque
from multiprocessing import Pool, Manager

from .._util import Capture
from .context import Context

calculator = None


def worker(calc_proxy, cxt):
    global calculator
    if calculator is None:
        calculator = calc_proxy[0]

    with Capture() as capture:
        r = list(calculator._calculate(cxt))

        return r, capture.result


class MolPool(object):
    def __init__(self, calc, nproc):
        self.pool = Pool(nproc)
        self.mgr = Manager()
        self.calc = calc
        self.calc_proxy = self.mgr.list([calc])
        self.nproc = nproc

    def __enter__(self):
        self.mgr.__enter__()
        return self

    def __exit__(self, *args, **kwargs):
        self.pool.terminate()
        self.mgr.__exit__(*args, **kwargs)

    def map(self, mols, id):
        return MolIterator(self, mols, id, self.nproc * 2 + 10)

    def submit(self, mol, id):
        cxt = Context.from_calculator(self.calc, mol, id)
        return self.pool.apply_async(worker, (self.calc_proxy, cxt))


class MolIterator(object):
    def __init__(self, pool, mols, id, buf):
        self.pool = pool
        self.futures = deque()
        self.mols = iter(mols)
        self.id = id

        for mol in islice(self.mols, buf):
            self.submit(mol)

    def submit(self, mol):
        self.futures.append((mol, self.pool.submit(mol, self.id)))

    def __iter__(self):
        return self

    def __next__(self):
        try:
            self.submit(next(self.mols))
        except StopIteration:
            pass

        try:
            mol, fut = self.futures.popleft()
            return mol, fut.get()
        except IndexError:
            raise StopIteration

    next = __next__


def parallel(calc, mols, nproc, nmols, quiet, ipynb, id):
    with MolPool(calc, nproc) as pool, calc._progress(quiet, nmols, ipynb) as bar:
        for mol, (r, err) in pool.map(mols, id):
            for e in err:
                e = e.rstrip()
                if not e:
                    continue

                bar.write(e)

            yield calc._wrap_result(mol, r)
            bar.update()
