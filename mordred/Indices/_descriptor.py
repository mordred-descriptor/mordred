from .._base import Descriptor
from ._eccentric_connectivity import EccentricConnectivityIndex
from ._zagreb import Zagreb
from ._wiener import Wiener
from ._detour import Detour
from ._topological import Radius, Diameter, TopologicalShapeIndex


_descriptors = [
    EccentricConnectivityIndex,
    Zagreb,
    Wiener,
    Detour,
    Radius, Diameter, TopologicalShapeIndex,
]
__all__ = [d.__name__ for d in _descriptors]
