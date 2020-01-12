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

#####matrix
```
Overlap.matrix(self, basis) -> np.array
```


Return a matrix of the operator projected onto a basis.





####Properties

#####hermitian = True





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

#####matrix
```
Kinetic.matrix(self, basis) -> np.array
```


Return a matrix of the operator projected onto a basis.





####Properties

#####hermitian = True





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

#####matrix
```
Hamiltonian.matrix(self, basis) -> np.array
```


Return a matrix of the operator projected onto a basis.





####Properties

#####hermitian = True





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

#####matrix
```
Potential.matrix(self, basis) -> np.array
```


Return a matrix of the operator projected onto a basis.





####Properties

#####hermitian = True





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

#####matrix
```
ExtendedHuckelHamiltonian.matrix(self, basis=None) -> np.array
```


Create the Hamiltonian under the Extended Hueckel Approximation.





####Properties

#####hermitian = True





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

