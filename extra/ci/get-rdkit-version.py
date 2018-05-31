def main(data, os, pyver):
    d = {
        (os, pyver): ver
        for os, pyver, ver
        in (line.strip().split() for line in open(data))
    }

    print(d[(os, pyver)])


if __name__ == "__main__":
    import sys
    import os

    data = os.path.join(os.path.dirname(__file__), "rdkit-versions.txt")
    main(data, *sys.argv[1:])
