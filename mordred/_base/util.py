from ..error import MissingValueBase


def is_missing(v):
    """Check argument is either MissingValue or not.

    >>> from mordred.error import Missing, Error
    >>> is_missing(1)
    False
    >>> is_missing(Missing(ValueError("missing"), ()))
    True
    >>> is_missing(Error(ValueError("error"), ()))
    True

    Parameters:
        v(any): value

    Returns:
        bool

    """
    return isinstance(v, MissingValueBase)
