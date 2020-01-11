# quantized.utils Module



## Classes

---

###Parabola

```
Parabola(a: 'float', b: 'float', c: 'float')
```
#### Fields

 **a (float)** 

> Constraints:  <function not_nan at 0x7fb228805620>, <function not_inf at 0x7fb228805598>, Range(min=-1000.0, max=1000.0)

 **b (float)** 

> Constraints:  <function not_nan at 0x7fb228805620>, <function not_inf at 0x7fb228805598>, Range(min=-1000.0, max=1000.0)

 **c (float)** 

> Constraints:  <function not_nan at 0x7fb228805620>, <function not_inf at 0x7fb228805598>, Range(min=-1000.0, max=1000.0)



####Properties

#####has\_vertex


None
#####vertex


None



####Static Methods

#####from\_points
```
Parabola.from_points(
    point1: 'Tuple[float, float]',
    point2: 'Tuple[float, float]',
    point3: 'Tuple[float, float]'
) -> Parabola
```


Create a Parabola passing through 3 points.

Parameters
-----------
point1
    x, y point
point2
    x, y point
point3
    x, y point

Returns
--------
Parabola
    Parabola with coefficients fit to the points.



####Dunder Methods

#####\_\_call\_\_
```
Parabola.__call__(self, x) -> <class 'inspect._empty'>
```


Call self as a function.



## Functions

----

###pairwise_array_from_func

```
pairwise_array_from_func(
    items: 'Sequence[T]',
    func: 'Callable[[T, T], float]',
    symmetric=False
) -> np.ndarray
```


Create a pairwise array by applying a function to all pairs of items.


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