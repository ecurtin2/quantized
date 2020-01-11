# first line: 354
@cache
def get_expansion_coeffs(state: Callable, basis: List[Callable]) -> List[float]:
    return [Overlap()(state, basis_func) for basis_func in basis]
