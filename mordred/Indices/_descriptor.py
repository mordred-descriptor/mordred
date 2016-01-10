from .._base import Descriptor
from ._eccentric_connectivity import EccentricConnectivityIndex
from ._zagreb import ZagrebIndex
from ._wiener import WienerIndex
from ._detour import DetourIndex
from ._topological import Radius, Diameter, TopologicalShapeIndex, PetitjeanIndex
from ._kier import KappaShapeIndex


_descriptors = [
    EccentricConnectivityIndex,
    ZagrebIndex,
    WienerIndex,
    DetourIndex,
    Radius, Diameter, TopologicalShapeIndex, PetitjeanIndex,
    KappaShapeIndex,
]
__all__ = [d.__name__ for d in _descriptors]
