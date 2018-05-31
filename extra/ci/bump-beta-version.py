from distutils.version import StrictVersion
import sys


def main(version):
    (major, minor, patch) = StrictVersion(version).version
    print("{}.{}.{}a1".format(major, minor, patch + 1))


if __name__ == "__main__":
    main(sys.argv[1])
