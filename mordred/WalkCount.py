from ._walk_count import WalkCount

__all__ = ('WalkCount',)

if __name__ == '__main__':
    from .__main__ import submodule
    submodule([
        WalkCount,
    ])
