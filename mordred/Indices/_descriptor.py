from .._base import Descriptor
from ._eccentric_connectivity import EccentricConnectivityIndex


_descriptors = [EccentricConnectivityIndex]
__all__ = [d.__name__ for d in _descriptors]
