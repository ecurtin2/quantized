# quantized.potentials Module



## Classes

---

###TripleWell

```
TripleWell(
    center1: Tuple[float, float],
    barrier12: Tuple[float, float],
    center2: Tuple[float, float],
    barrier23: Tuple[float, float],
    center3: Tuple[float, float],
    well1: quantized.utils.Parabola,
    well2: quantized.utils.Parabola,
    well3: quantized.utils.Parabola
)
```
Construct a well potential from the points.
```
x                                                         x
x                                                         x
x                                                         x
x                                                        x
x                                                        x
 x          barrier12                                   x
 x         x         x            barrier23            x
  x       x          x           x         x          x
   center1            x         x           x       x
                       x       x             center3
                        center2
```



#### Fields

 **center1 (float | float)** 

 **barrier12 (float | float)** 

 **center2 (float | float)** 

 **barrier23 (float | float)** 

 **center3 (float | float)** 

 **well1 (Parabola)** 

 **well2 (Parabola)** 

 **well3 (Parabola)** 







####Static Methods

#####from\_params
```
TripleWell.from_params(
    well1_depth: float,
    well1_halfwidth: float,
    bridge_length: float,
    bridge_depth: float,
    well3_halfwidth: float,
    well3_depth: float
) -> <class 'inspect._empty'>
```






####Dunder Methods

#####\_\_call\_\_
```
TripleWell.__call__(
    self,
    x: Union[numpy.ndarray, float]
) -> typing.Union[numpy.ndarray, float]
```


Call self as a function.

 --- 

###Harmonic

```
Harmonic(center: float, mass: float = 1.0, omega: float = 1.0)
```
#### Fields

 **center (float)** 

> Constraints:  Range(min=-1000.0, max=1000.0)

 **mass (float)** 

> Default: 1.0

> Constraints:  Range(min=1e-08, max=1000.0)

 **omega (float)** 

> Default: 1.0

> Constraints:  Range(min=1e-08, max=1000.0)









####Dunder Methods

#####\_\_call\_\_
```
Harmonic.__call__(self, x: float) -> <class 'float'>
```


Return the value of the harmonic potential at coordinate x

