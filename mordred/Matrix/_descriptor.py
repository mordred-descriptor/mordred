from ._barysz import BaryszMatrix
from ._burden import BCUT
from ._detour import DetourMatrix


_descriptors = [BaryszMatrix, BCUT, DetourMatrix]
__all__ = [d.__name__ for d in _descriptors]
