from .._base import Descriptor
from ._eccentric_connectivity import EccentricConnectivityIndex
from ._zagreb import Zagreb


_descriptors = [EccentricConnectivityIndex, Zagreb]
__all__ = [d.__name__ for d in _descriptors]
