from ._barysz import BaryszMatrix
from ._burden import BCUT
from ._detour import DetourMatrix
from ._distance import DistanceMatrix


_descriptors = [
    BaryszMatrix,
    BCUT,
    DetourMatrix,
    DistanceMatrix,
]

__all__ = [d.__name__ for d in _descriptors]
