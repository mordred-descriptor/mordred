from ._information_content import (
    InformationContent,
    TotalIC, StructuralIC, BondingIC, ComplementaryIC,
    ModifiedIC, ZModifiedIC,
)

__all__ = (
    'InformationContent',
    'TotalIC', 'StructuralIC', 'BondingIC', 'ComplementaryIC',
    'ModifiedIC', 'ZModifiedIC',
)

if __name__ == '__main__':
    from .__main__ import submodule
    submodule()
