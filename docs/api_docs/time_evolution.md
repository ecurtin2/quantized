# quantized.time_evolution Module



## Classes

---

###TimeEvolvingState

```
TimeEvolvingState(
    initial_state: Callable,
    eigen_basis: quantized.basis.EigenBasis
)
```
#### Fields

 **initial_state ()** 

 **eigen_basis (EigenBasis)** 

####Methods

#####observable
```
TimeEvolvingState.observable(
    self,
    operator: Callable,
    hermitian: bool
) -> TimeEvolvingObservable
```












####Dunder Methods

#####\_\_call\_\_
```
TimeEvolvingState.__call__(
    self,
    x: Union[float, numpy.ndarray],
    t: float
) -> <class 'float'>
```


Call self as a function.

 --- 

###TimeEvolvingObservable

```
TimeEvolvingObservable(
    time_evolving_state: quantized.time_evolution.TimeEvolvingState,
    operator: quantized.operators.Operator,
    hermitian: bool
)
```
#### Fields

 **time_evolving_state (TimeEvolvingState)** 

 **operator (Operator)** 

 **hermitian (bool)** 









####Dunder Methods

#####\_\_call\_\_
```
TimeEvolvingObservable.__call__(
    self,
    t: Union[float, <built-in function array>]
) -> typing.Union[float, <built-in function array>]
```


Call self as a function.

