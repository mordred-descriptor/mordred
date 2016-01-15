from ._matrix import DetourMatrix
from ._index import DetourIndex


_descriptors = [
    DetourMatrix,
    DetourIndex,
]

__all__ = [d.__name__ for d in _descriptors]
