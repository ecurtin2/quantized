def check_range_inclusive(name: str, x, min_val=None, max_val=None) -> str:
    """Check if value is between min and max values, returning an error message if not.

    min_val and max_val must be comparable to x


    Parameters
    -----------
    name
        Name to give the variable in the error message
    x
        The value
    min_val
        Optional minimum value. Inclusive
    max_val
        Optional max value. Inclusive


    Returns
    --------
    error
        A string containing an error message if the variable
        is out of range. If there is no error, it will be an
        empty string.

    Examples
    ---------

    >>> from transit_chem.validation import check_range_inclusive
    >>> check_range_inclusive("x", 3, 0, 2)
    'x out of range, got 3, expected 0 < x < 2'
    >>> check_range_inclusive("x", 3, 0, 10)
    ''
    >>> check_range_inclusive("x", 3, -1)
    ''
    >>> check_range_inclusive("x", -1, 2)
    'x too small, got -1, expected x > 2'
    >>> check_range_inclusive("x", 10, max_val=10)
    ''
    >>> check_range_inclusive("x", 15, max_val=10)
    'x too large, got 15, expected x < 10'

    """
    if (min_val is None) and (max_val is None):
        raise ValueError("Min and max_val can't both be None")
    elif min_val is None:
        if x > max_val:
            return f"{name} too large, got {x}, expected {name} < {max_val}"
    elif max_val is None:
        if x < min_val:
            return f"{name} too small, got {x}, expected {name} > {min_val}"
    else:
        if x < min_val or x > max_val:
            return (
                f"{name} out of range, got {x}, expected {min_val} < {name} < {max_val}"
            )

    return ""
