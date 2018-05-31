from __future__ import print_function
import sys
import re


def main(path):
    contents = [line.rstrip() for line in open(path)]
    with open(path, 'w') as dst:
        for line in contents:
            if '-' == line[0]:
                continue

            vs = re.split(r'\s*;\s*', line)
            if len(vs) == 1:
                print(line, file=dst)
                continue

            if len(vs) > 2:
                raise ValueError("#field > 2")

            v = sys.version_info
            pkg, cond = vs
            cond = cond.replace("python_version", '{}.{}'.format(v.major, v.minor))
            cond = re.sub(r"'(\d+\.\d+)'", r'\1', cond)
            if eval(cond):
                print(pkg, file=dst)


if __name__ == "__main__":
    main(sys.argv[1])
