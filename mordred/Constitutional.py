from ._constitutional import ConstitutionalSum, ConstitutionalMean

__all__ = 'ConstitutionalSum', 'ConstitutionalMean',

if __name__ == '__main__':
    from .__main__ import submodule
    submodule([
        ConstitutionalSum, ConstitutionalMean
    ])
