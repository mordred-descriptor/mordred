from ._detour_matrix import DetourMatrix
from ._detour_index import DetourIndex

__all__ = 'DetourMatrix', 'DetourIndex',

if __name__ == '__main__':
    from .__main__ import submodule
    submodule([
        DetourMatrix,
        DetourIndex,
    ])
