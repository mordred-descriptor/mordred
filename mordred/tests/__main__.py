import os

import nose


def main():
    base = os.path.dirname(os.path.dirname(__file__))
    hidden = [
        os.path.join(base, n)
        for n in os.listdir(base)
        if n[:1] == "_" and os.path.splitext(n)[1] == ".py"
    ]

    tests = [base, os.path.join(base, "_base")] + hidden

    os.environ["NOSE_WITH_DOCTEST"] = "1"

    nose.main(defaultTest=",".join(tests))


if __name__ == "__main__":
    main()
