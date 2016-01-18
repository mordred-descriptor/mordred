from ._autocorrelation import AATS, AATSC, ATS, ATSC, GATS, MATS

__all__ = ('ATS', 'AATS', 'ATSC', 'AATSC', 'MATS', 'GATS',)

if __name__ == '__main__':
    from .__main__ import submodule
    submodule([
        ATS, AATS, ATSC, AATSC, MATS, GATS,
    ])
