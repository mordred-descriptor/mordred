from ._adjacency_matrix import AdjacencyMatrix

__all__ = 'AdjacencyMatrix',

if __name__ == '__main__':
    from .__main__ import submodule
    submodule([
        AdjacencyMatrix,
    ])
