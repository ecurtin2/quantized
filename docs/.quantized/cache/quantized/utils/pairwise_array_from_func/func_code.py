# first line: 48
@cache
def pairwise_array_from_func(
    items: Sequence[T], func: Callable[[T, T], float], symmetric=False
) -> np.ndarray:
    """Create a pairwise array by applying a function to all pairs of items.


    Parameters
    -----------
    items
        A container from which pairs will be generated. Must support len() and integer indexing over range(len(items))
    func
        A function f(first, second, *args, **kwargs) which takes 2 items and returns a float.
    symmetric
        Whether the resulting matrix should be symmetric. If true,
        will only compute each (i, j) pair once and set both [i, j] and [j, i] to that value.

    Returns
    --------
    np.array
        The resulting matrix

    Examples
    ---------

    >>> from quantized.utils import pairwise_array_from_func
    >>> def distance(i, j):
    ...     return abs(i - j)
    ...
    >>> pairwise_array_from_func([1, 2, 4], distance)
    array([[0., 1., 3.],
           [1., 0., 2.],
           [3., 2., 0.]])
    >>> pairwise_array_from_func([1, 2, 4], distance, symmetric=True)
    array([[0., 1., 3.],
           [1., 0., 2.],
           [3., 2., 0.]])


    """
    n = len(items)
    result = np.zeros((n, n))

    fut_res: Dict[Future, Tuple] = {}

    if symmetric:
        combs = list(combinations_with_replacement(range(n), 2))
    else:
        combs = list(product(range(n), repeat=2))

    with ProcessPoolExecutor() as exc:
        for i, j in combs:
            fut_res[exc.submit(func, items[i], items[j])] = (i, j)

        with tqdm(total=len(combs), disable=not conf.enable_progressbar, ncols=100) as pbar:
            for future in as_completed(fut_res.keys()):
                res = future.result()
                i, j = fut_res[future]
                if symmetric:
                    result[i, j] = result[j, i] = res
                else:
                    result[i, j] = res
                pbar.update(1)

    return result
