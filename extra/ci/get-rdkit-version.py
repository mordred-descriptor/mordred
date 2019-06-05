import re

comment = re.compile(r'^\s*#')


def main(data, os, pyver):
    d = {
        (os, pyver): ver
        for os, pyver, ver
        in (line.strip().split() for line in open(data) if not comment.match(line))
    }

    print(d[(os, pyver)])


if __name__ == "__main__":
    import sys

    main(*sys.argv[1:])
