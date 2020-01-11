# quantized.hopping_matrix Module



## Classes

---

###OccupancyProbabilites

```
OccupancyProbabilites(s)
```
#### Fields

 **s (TimeEvolvingObservable)** 



####Properties

#####initial


None



####Static Methods

#####from\_1d\_state
```
OccupancyProbabilites.from_1d_state(
    state: quantized.time_evolution.TimeEvolvingState,
    borders: Tuple[float, ...]
) -> OccupancyProbabilites
```


None



####Dunder Methods

#####\_\_call\_\_
```
OccupancyProbabilites.__call__(self, t: float) -> <built-in function array>
```


Call self as a function.
#####\_\_getitem\_\_
```
OccupancyProbabilites.__getitem__(self, item) -> typing.Callable
```


None
#####\_\_iter\_\_
```
OccupancyProbabilites.__iter__(self) -> <class 'inspect._empty'>
```


None

 --- 

###HoppingMatrix

```
HoppingMatrix(occ_probs: quantized.hopping_matrix.OccupancyProbabilites)
```
#### Fields

 **occ_probs (OccupancyProbabilites)** 

####Methods

#####N


int([x]) -> integer
int(x, base=10) -> integer

Convert a number or string to an integer, or return 0 if no arguments
are given.  If x is a number, return x.__int__().  For floating point
numbers, this truncates towards zero.

If x is not a number or if base is given, then x must be a string,
bytes, or bytearray instance representing an integer literal in the
given base.  The literal can be preceded by '+' or '-' and be surrounded
by whitespace.  The base defaults to 10.  Valid bases are 0 and 2-36.
Base 0 means to interpret the base from the string as an integer literal.
>>> int('0b100', base=0)
4
#####at\_time
```
HoppingMatrix.at_time(self, t: float, delta_t: float) -> <built-in function array>
```


None







####Dunder Methods

#####\_\_call\_\_
```
HoppingMatrix.__call__(
    self,
    t: Union[float, numpy.ndarray],
    delta_t: float
) -> <class 'inspect._empty'>
```


Call self as a function.

 --- 

###Pnot

```
Pnot(
    acceptor: int,
    hopping_matrix: quantized.hopping_matrix.HoppingMatrix,
    times: List[float],
    values: List[float],
    delta_t: float
)
```
#### Fields

 **acceptor (int)** 

 **hopping_matrix (HoppingMatrix)** 

 **times (float)** 

 **values (float)** 

 **delta_t (float)** 

####Methods

#####time\_when\_equal\_to


None



####Properties

#####acceptor\_occ\_prob


None
#####final\_value


None
#####initial\_value


None
#####max\_time


None
#####min\_time


None
#####tau90


None



####Static Methods

#####converged\_with\_timestep
```
Pnot.converged_with_timestep(
    a: quantized.hopping_matrix.HoppingMatrix,
    acceptor: int,
    until_equals: float,
    tolerance: float,
    max_dt: float = 1.0,
    min_dt: float = 0.0001
) -> Pnot
```


None
#####gen
```
Pnot.gen(
    delta_t: float,
    hopping_matrix: quantized.hopping_matrix.HoppingMatrix,
    acceptor: int
) -> typing.Generator[typing.Tuple, NoneType, NoneType]
```


None
#####gen\_until\_prob
```
Pnot.gen_until_prob(
    p: float,
    delta_t: float,
    hopping_matrix: quantized.hopping_matrix.HoppingMatrix,
    acceptor: int
) -> typing.Iterator[typing.Tuple[float, float]]
```


None
#####gen\_until\_time
```
Pnot.gen_until_time(
    t: float,
    delta_t: float,
    hopping_matrix: quantized.hopping_matrix.HoppingMatrix,
    acceptor: int
) -> typing.Iterator[typing.Tuple[float, float]]
```


None



####Dunder Methods

#####\_\_call\_\_
```
Pnot.__call__(self, t: float) -> <class 'float'>
```


Call self as a function.



## Functions

----

