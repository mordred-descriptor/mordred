from ._topological_index import Diameter, Radius, TopologicalShapeIndex, PetitjeanIndex

__all__ = 'Diameter', 'Radius', 'TopologicalShapeIndex', 'PetitjeanIndex',

if __name__ == '__main__':
    from .__main__ import submodule
    submodule([
        Diameter, Radius, TopologicalShapeIndex, PetitjeanIndex,
    ])
