from ._autocorrelation import ATS, AATS, ATSC, AATSC, MATS, GATS

if __name__ == '__main__':
    from .__main__ import submodule
    submodule([
        ATS, AATS, ATSC, AATSC, MATS, GATS,
    ])
