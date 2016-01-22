r"""Wildman-Crippen LogP/MR descriptor.

References
    * :cite:`10.1021/ci990307l`
"""

from ._wildman_crippen_logp import WildmanCrippenLogP, WildmanCrippenMR

__all__ = (
    'WildmanCrippenLogP',
    'WildmanCrippenMR',
)

if __name__ == '__main__':
    from .__main__ import submodule
    submodule()
