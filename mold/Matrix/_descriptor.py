from ._barysz import BaryszMatrix
from ._burden import BCUT


_descriptors = [BaryszMatrix, BCUT]
__all__ = [d.__name__ for d in _descriptors]
