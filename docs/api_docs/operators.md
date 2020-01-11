# quantized.operators Module



## Classes

---

###Overlap

```
Overlap(lower_limit: 'float' = -inf, upper_limit: 'float' = inf)
```
#### Fields

 **lower_limit (float)** 

> Default: -inf

 **upper_limit (float)** 

> Default: inf

####Methods

#####hermitian


bool(x) -> bool

Returns True when the argument x is true, False otherwise.
The builtins True and False are the only two instances of the class bool.
The class bool is a subclass of the class int, and cannot be subclassed.
#####matrix
```
Overlap.matrix(self, basis) -> np.array
```


None







####Dunder Methods

#####\_\_call\_\_
```
Overlap.__call__(self, first: 'Callable', second: 'Callable') -> float
```


Call self as a function.

 --- 

###Kinetic

```
Kinetic(
    overlap: 'Overlap' = Overlap(lower_limit=-inf, upper_limit=inf)
)
```
#### Fields

 **overlap (Overlap)** 

> Default: Overlap(lower_limit=-inf, upper_limit=inf)

####Methods

#####hermitian


bool(x) -> bool

Returns True when the argument x is true, False otherwise.
The builtins True and False are the only two instances of the class bool.
The class bool is a subclass of the class int, and cannot be subclassed.
#####matrix
```
Kinetic.matrix(self, basis) -> np.array
```


None







####Dunder Methods

#####\_\_call\_\_
```
Kinetic.__call__(self, first, second) -> float
```


Call self as a function.

 --- 

###Hamiltonian

```
Hamiltonian(
    potential: 'Callable[[float], float]',
    kinetic: 'Kinetic' = Kinetic(overlap=Overlap(lower_limit=-inf, upper_limit=inf))
)
```
#### Fields

 **potential (Callable[[float], float])** 

 **kinetic (Kinetic)** 

> Default: Kinetic(overlap=Overlap(lower_limit=-inf, upper_limit=inf))

####Methods

#####hermitian


bool(x) -> bool

Returns True when the argument x is true, False otherwise.
The builtins True and False are the only two instances of the class bool.
The class bool is a subclass of the class int, and cannot be subclassed.
#####matrix
```
Hamiltonian.matrix(self, basis) -> np.array
```


None







####Dunder Methods

#####\_\_call\_\_
```
Hamiltonian.__call__(self, first, second) -> float
```


Call self as a function.

 --- 

###Potential

```
Potential(
    potential: 'Callable[[float], float]',
    overlap: 'Overlap' = Overlap(lower_limit=-inf, upper_limit=inf)
)
```
#### Fields

 **potential (Callable[[float], float])** 

 **overlap (Overlap)** 

> Default: Overlap(lower_limit=-inf, upper_limit=inf)

####Methods

#####hermitian


bool(x) -> bool

Returns True when the argument x is true, False otherwise.
The builtins True and False are the only two instances of the class bool.
The class bool is a subclass of the class int, and cannot be subclassed.
#####matrix
```
Potential.matrix(self, basis) -> np.array
```


None







####Dunder Methods

#####\_\_call\_\_
```
Potential.__call__(self, first, second) -> float
```


Call self as a function.

 --- 

###ExtendedHuckelHamiltonian

```
ExtendedHuckelHamiltonian(S: 'np.array', molecule: 'Molecule')
```
#### Fields

 **S (np.array)** 

 **molecule (Molecule)** 

####Methods

#####hermitian


bool(x) -> bool

Returns True when the argument x is true, False otherwise.
The builtins True and False are the only two instances of the class bool.
The class bool is a subclass of the class int, and cannot be subclassed.
#####matrix
```
ExtendedHuckelHamiltonian.matrix(self, basis=None) -> np.array
```


Create the Hamiltonian under the Extended Hueckel Approximation.







####Dunder Methods

#####\_\_call\_\_
```
ExtendedHuckelHamiltonian.__call__(
    self,
    first: 'Callable[[float], float]',
    second: 'Callable[[float], float]'
) -> float
```


Call self as a function.



## Functions

----

