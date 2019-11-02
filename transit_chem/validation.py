def check_range_inclusive(name: str, x, min_val=None, max_val=None) -> str:
    
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
                f"{name} out of range, got {x}" f"expected {min_val} < {name} < {max_val}"
            )

    return ""
