from .._base import Descriptor
from ._eccentric_connectivity import EccentricConnectivityIndex
from ._zagreb import Zagreb
from ._wiener import Wiener


_descriptors = [EccentricConnectivityIndex, Zagreb, Wiener]
__all__ = [d.__name__ for d in _descriptors]
