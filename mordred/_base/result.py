from .util import is_missing


class Result(list):
    r"""Result type."""

    def __init__(self, r, d):
        super(Result, self).__init__(r)
        self._descriptors = d

    def fillna(self, value=float("nan")):
        r"""Replace missing value to "value".

        Parameters:
            value: value that missing value is replaced

        Returns:
            Result

        """
        return self.__class__(
            [(value if is_missing(v) else v) for v in self],
            self._descriptors,
        )

    def dropna(self):
        r"""Delete missing value.

        Returns:
            Result

        """
        newvalues = []
        newdescs = []
        for v, d in zip(self, self._descriptors):
            if not is_missing(v):
                newvalues.append(v)
                newdescs.append(d)

        return self.__class__(newvalues, newdescs)

    def asdict(self, rawkey=False):
        r"""Convert Result to dict.

        Parameters:
            rawkey(bool):
                * True: dict key is Descriptor instance
                * False: dict key is str

        Returns:
            dict

        """
        if rawkey:
            def keyconv(k):
                return k
        else:
            keyconv = str

        return {
            keyconv(k): v
            for k, v in zip(self._descriptors, self)
        }
