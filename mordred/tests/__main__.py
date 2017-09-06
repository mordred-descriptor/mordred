import os

import nose


def main():
    base = os.path.dirname(os.path.dirname(__file__))
    tests = [base, os.path.join(base, "_base")]

    os.environ["NOSE_WITH_DOCTEST"] = "1"

    nose.main(
        defaultTest=",".join(tests),
    )


if __name__ == "__main__":
    main()
