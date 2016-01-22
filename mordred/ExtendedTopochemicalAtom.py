r"""Extended Topochemical Atom(ETA) descriptor.

References
    * :cite:`10.1021/ci0342066`
    * :cite:`10.1016/j.jhazmat.2013.03.023`
"""

from ._extended_topochemical_atom import (
    EtaCoreCount, EtaShapeIndex,
    EtaVEMCount,
    EtaCompositeIndex, EtaFunctionalityIndex, EtaBranchingIndex,

    EtaDeltaAlpha,
    EtaEpsilon, EtaDeltaEpsilon,
    EtaDeltaBeta,
    EtaPsi, EtaDeltaPsi,
)

__all__ = (
    'EtaCoreCount', 'EtaShapeIndex',
    'EtaVEMCount',
    'EtaCompositeIndex', 'EtaFunctionalityIndex', 'EtaBranchingIndex',

    'EtaDeltaAlpha',
    'EtaEpsilon', 'EtaDeltaEpsilon',
    'EtaDeltaBeta',
    'EtaPsi', 'EtaDeltaPsi',
)

if __name__ == '__main__':
    from .__main__ import submodule

    submodule([
        EtaCoreCount, EtaShapeIndex,
        EtaVEMCount,
        EtaCompositeIndex, EtaFunctionalityIndex, EtaBranchingIndex,

        EtaDeltaAlpha,
        EtaEpsilon, EtaDeltaEpsilon,
        EtaDeltaBeta,
        EtaPsi, EtaDeltaPsi,
    ])
