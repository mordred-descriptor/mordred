r"""Wildman-Crippen LogP/MR descriptor.

References
    * :cite:`10.1021/ci990307l`
"""

from ._slogp import SLogP, SMR

__all__ = ('SLogP', 'SMR',)

if __name__ == '__main__':
    from .__main__ import submodule
    submodule()
