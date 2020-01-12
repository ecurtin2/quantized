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








####Static Methods

#####from\_1d\_state
```
OccupancyProbabilites.from_1d_state(
    state: quantized.time_evolution.TimeEvolvingState,
    borders: Tuple[float, ...]
) -> OccupancyProbabilites
```






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



#####\_\_iter\_\_
```
OccupancyProbabilites.__iter__(self) -> <class 'inspect._empty'>
```




 --- 

###HoppingMatrix

```
HoppingMatrix(occ_probs: quantized.hopping_matrix.OccupancyProbabilites)
```
#### Fields

 **occ_probs (OccupancyProbabilites)** 

####Methods

#####at\_time
```
HoppingMatrix.at_time(self, t: float, delta_t: float) -> <built-in function array>
```








####Properties

#####N = 3





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






####Properties

#####acceptor\_occ\_prob



#####final\_value



#####initial\_value



#####max\_time



#####min\_time



#####tau90








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



#####gen
```
Pnot.gen(
    delta_t: float,
    hopping_matrix: quantized.hopping_matrix.HoppingMatrix,
    acceptor: int
) -> typing.Generator[typing.Tuple, NoneType, NoneType]
```



#####gen\_until\_prob
```
Pnot.gen_until_prob(
    p: float,
    delta_t: float,
    hopping_matrix: quantized.hopping_matrix.HoppingMatrix,
    acceptor: int
) -> typing.Iterator[typing.Tuple[float, float]]
```



#####gen\_until\_time
```
Pnot.gen_until_time(
    t: float,
    delta_t: float,
    hopping_matrix: quantized.hopping_matrix.HoppingMatrix,
    acceptor: int
) -> typing.Iterator[typing.Tuple[float, float]]
```






####Dunder Methods

#####\_\_call\_\_
```
Pnot.__call__(self, t: float) -> <class 'float'>
```


Call self as a function.

