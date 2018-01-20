import numpy as np
from six import string_types, integer_types

from .util import is_missing
from .descriptor import Descriptor


class Result(object):
    r"""Result type."""

    __slots__ = ("_values", "_descriptors", "_name_to_index")

    def __init__(self, r, d):
        self._values = list(r)
        self._descriptors = list(d)
        self._name_to_index = {str(a): i for i, a in enumerate(d)}

    def __str__(self):
        return str(self._name_to_index)

    def __repr__(self):
        return "{}({!r},{!r})".format(
            self.__class__.__name__,
            self._values,
            self._descriptors,
        )

    def fill_missing(self, value=np.nan):
        r"""Replace missing value to "value".

        Parameters:
            value: value that missing value is replaced

        Returns:
            Result

        """
        return self.__class__(
            [(value if is_missing(v) else v) for v in self.values()],
            self.keys(),
        )

    def drop_missing(self):
        r"""Delete missing value.

        Returns:
            Result

        """
        newvalues = []
        newdescs = []
        for d, v in self.items():
            if not is_missing(v):
                newvalues.append(v)
                newdescs.append(d)

        return self.__class__(newvalues, newdescs)

    def items(self):
        r"""Get items.

        Returns:
            Iterable[(Descriptor, value)]

        """
        return ((k, v) for k, v in zip(self.keys(), self.values()))

    def keys(self):
        r"""Get descriptors instances.

        Returns:
            Iterable[Descriptor]

        """
        return iter(self._descriptors)

    def values(self):
        r"""Get descriptor values.

        Returns:
            Iterable[value]

        """
        return iter(self._values)

    __iter__ = values

    def __reversed__(self):
        return reversed(self._values)

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
            return dict(self.items())
        else:
            return {
                str(k): v
                for k, v in self.items()
            }

    @property
    def ix(self):
        r"""Access descriptor value by index.

        >>> from mordred import Calculator, Lipinski
        >>> from rdkit import Chem
        >>> result = Calculator(Lipinski.Lipinski)(Chem.MolFromSmiles("C1CCCCC1"))
        >>> result.ix[0]
        True
        """
        return GetValueByIndex(self._values)

    @property
    def name(self):
        r"""Access descriptor value by descriptor name or instance.

        >>> from mordred import Calculator, descriptors
        >>> from rdkit import Chem
        >>> result = Calculator(descriptors)(Chem.MolFromSmiles("C1CCCCC1"))
        >>> result.name["C2SP3"]
        6

        """
        return GetValueByName(self._values, self._name_to_index)

    def __getitem__(self, key):
        if isinstance(key, (integer_types, slice)):
            return self.ix[key]

        elif isinstance(key, (string_types, Descriptor)):
            return self.name[key]

        else:
            raise TypeError(
                "Result indices must be "
                "integers, slices, strings or Descriptor instance, "
                "not {}".format(key.__class__.__name__))

    def __len__(self):
        return len(self._descriptors)


class GetValueByIndex(object):
    __slots__ = ("_values",)

    def __init__(self, values):
        self._values = values

    def __getitem__(self, key):
        return self._values[key]


class GetValueByName(object):
    __slots__ = ("_values", "_name_to_index")

    def __init__(self, values, name_to_index):
        self._values = values
        self._name_to_index = name_to_index

    def __getitem__(self, key):
        return self._values[self._name_to_index[str(key)]]
