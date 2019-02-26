import sys
from distutils.version import StrictVersion


def main(version):
    vsn = StrictVersion(version)
    if vsn.prerelease:
        print(version)
        return

    (major, minor, patch) = vsn.version
    print("{}.{}.{}a1".format(major, minor, patch + 1))


if __name__ == "__main__":
    main(sys.argv[1])
