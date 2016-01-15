from ._weight import Weight
from ._mc_gowan_volume import McGowanVolume
from ._vdw_volume_abc import VdwVolumeABC
from ._wildman_crippen_logp import WildmanCrippenLogP

_descriptors = [
    Weight,
    McGowanVolume,
    VdwVolumeABC,
    WildmanCrippenLogP,
]

__all__ = [d.__name__ for d in _descriptors]
