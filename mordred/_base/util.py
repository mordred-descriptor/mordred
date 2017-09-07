from ..error import MissingValueBase


def is_missing(v):
    """Check argument is either MissingValue or not.

    Parameters:
        v(any): value

    Returns:
        bool

    """
    return isinstance(v, MissingValueBase)
