r"""charged partial surface area descriptor.

References
    * :cite:`10.1021/ac00220a013`
"""

from ._cpsa import (
    PNSA, PPSA,
    DPSA,
    FNSA, FPSA,
    WNSA, WPSA,
    RNCG, RPCG,
    RNCS, RPCS,
    TASA, TPSA,
    RASA, RPSA,
)

__all__ = (
    'PNSA', 'PPSA',
    'DPSA',
    'FNSA', 'FPSA',
    'WNSA', 'WPSA',
    'RNCG', 'RPCG',
    'RNCS', 'RPCS',
    'TASA', 'TPSA',
    'RASA', 'RPSA',
)

if __name__ == '__main__':
    from .__main__ import submodule
    submodule()
