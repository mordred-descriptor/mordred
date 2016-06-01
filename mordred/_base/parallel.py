import six
from multiprocessing import Pool
from .context import Context


calculator = None


def initializer(self):
    global calculator

    calculator = self


def worker(cxt):
    return list(calculator._calculate(cxt))


def parallel(self, mols, nproc, callback):
    if six.PY3:
        def get_result(r):
            return r.get()
    else:
        # timeout: avoid python2 KeyboardInterrupt bug.
        # http://stackoverflow.com/a/1408476
        def get_result(r):
            return r.get(1e9)

    # without with-statement for compat. python2
    try:
        pool = Pool(nproc, initializer=initializer, initargs=(self,))

        def do_task(mol):
            args = Context.from_calculator(self, mol, -1)
            return pool.apply_async(worker, (args,))

        for m, result in [(m, do_task(m)) for m in mols]:
            yield m, get_result(result)
    finally:
        pool.terminate()
        pool.join()
