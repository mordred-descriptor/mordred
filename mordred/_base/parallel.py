from multiprocessing import Pool

import six

from .._util import Capture
from .context import Context

calculator = None


def initializer(self):
    global calculator

    calculator = self


def worker(cxt):
    with Capture() as capture:
        r = list(calculator._calculate(cxt))

        return r, capture.result


def parallel(self, mols, nproc, nmols, quiet, ipynb, id):
    if six.PY3:
        def get_result(r):
            return r.get()
    else:
        # timeout: avoid python2 KeyboardInterrupt bug.
        # http://stackoverflow.com/a/1408476
        def get_result(r):
            return r.get(1e9)

    # without with-statement for compat. python2
    pool = Pool(nproc, initializer=initializer, initargs=(self,))

    def do_task(mol):
        args = Context.from_calculator(self, mol, id)
        return pool.apply_async(worker, (args,))

    try:
        with self._progress(quiet, nmols, ipynb) as bar:
            for result in [do_task(m) for m in mols]:
                r, err = get_result(result)

                for e in err:
                    e = e.rstrip()
                    if not e:
                        continue

                    bar.write(e)

                yield self._wrap_result(r)
                bar.update()

    finally:
        pool.terminate()
        pool.join()
