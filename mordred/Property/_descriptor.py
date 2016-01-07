from ._weight import Weight
from ._mc_gowan_volume import McGowanVolume

_descriptors = [Weight, McGowanVolume]
__all__ = [d.__name__ for d in _descriptors]
