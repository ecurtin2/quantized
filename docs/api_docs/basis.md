# quantized.basis Module



## Classes

---

###HarmonicOscillator

```
HarmonicOscillator(
    n: 'int',
    center: 'float',
    mass: 'float' = 1.0,
    omega: 'float' = 1.0
)
```
A 1D quantum harmonic oscillator wave function.
```
>>> ho = HarmonicOscillator(n=2, center=0.5)
>>> ho
HarmonicOscillator(n=2, center=0.5, omega=1, mass=1.0)
>>> ho(0.5)
-0.5311259660135984
>>> ho(1000)
0.0
>>> ho(-1000)
0.0
```


#### Fields

 **n (int)** The quantum number

> Constraints:  Range(min=0, max=50)

 **center (float)** The center of the function

> Constraints:  Range(min=-1000.0, max=1000.0)

 **mass (float)** Mass of the particle

> Default: 1.0

> Constraints:  Range(min=1e-08, max=1000.0)

 **omega (float)** Angular frequency of the oscillator

> Default: 1.0

> Constraints:  Range(min=1e-08, max=1000.0)



####Properties

#####N


The normalization constant
#####energy


The energy of harmonic oscillator
#####potential


The potential for this oscillator





####Static Methods

#####from\_parabola
```
HarmonicOscillator.from_parabola(
    p: 'Parabola',
    n: 'int',
    mass: 'float' = 1.0
) -> HarmonicOscillator
```


Create a harmonic oscillator, who's potential is defined by the given parabola
#####from\_potential\_points
```
HarmonicOscillator.from_potential_points(
    point1: 'Tuple[float, float]',
    point2: 'Tuple[float, float]',
    point3: 'Tuple[float, float]',
    n: 'int',
    mass: 'float' = 1.0
) -> HarmonicOscillator
```


Create a harmonic oscillator wave function from 3 samples of the potential.

The three points are fit to a parabola (harmonic potential), then the parameters
for the harmonic oscillator are determined, and the corresponding wave function
generated and returned.

**point1**
A sample point of the potential

**point2**
A sample point of the potential

**point3**
A sample point of the potential

**n**
The quantum number of the resulting wave function

**mass**
The mass of the particle

**Examples**

```python
ho = HarmonicOscillator.from_potential_points(
...     point1=(0.5, 1),
...     point2=(2.0, 0.5),
...     point3=(3.0, 1.5),
...     n=0
... )
ho
HarmonicOscillator(n=0, center=1.5624999999999998, omega=1.0327955589886444, mass=1.0)
```



####Dunder Methods

#####\_\_call\_\_
```
HarmonicOscillator.__call__(
    self,
    x: 'Union[float, np.ndarray]'
) -> Union[float, np.ndarray]
```


Return
#####\_\_kinetic\_\_
```
HarmonicOscillator.__kinetic__(self) -> Callable
```


Return kinetic energy operator applied on this.
#####\_\_overlap\_\_
```
HarmonicOscillator.__overlap__(
    self,
    other,
    lower_limit: 'float',
    upper_limit: 'float'
) -> float
```


Determine the overlap with some other function

This specializes a generic overlap integral, and short circuits integral
calculations if the integral is analytically known.

 --- 

###EigenBasis

```
EigenBasis(
    states,
    energies,
    ao_S: 'np.ndarray',
    eigvecs: 'np.ndarray',
    ao_basis: 'List[Callable]'
)
```
A class for representing an eigenbasis for a hamiltonian

#### Fields

 **states (Tuple[Callable, ...])** A set of eigen states

 **energies (Tuple[float, ...])** The energies of the eigen states

 **ao_S (np.ndarray)** The overlap matrix in the original basis

 **eigvecs (np.ndarray)** The eigenvectors of the hamiltonian. Each column is a vector.

 **ao_basis (List[Callable])** The original basis

####Methods

#####transformed
```
EigenBasis.transformed(self, matrix: 'np.ndarray') -> np.ndarray
```


Given a matrix in the original basis, return the matrix in the Eigen basis.







####Static Methods

#####from\_basis
```
EigenBasis.from_basis(
    basis: 'List[Callable]',
    H: 'np.ndarray',
    S: 'np.ndarray'
) -> EigenBasis
```


Create an eigenbasis from another basis, given a hamiltonian and overlap matrix

**H (np.ndarray)**
The hamiltonian matrix in the basis

**S (np.ndarray)**
The overlap matrix in the basis



####Dunder Methods

#####\_\_len\_\_
```
EigenBasis.__len__(self) -> <class 'inspect._empty'>
```


The size of the eigenbasis



## Functions

----

###harmonic_basis_from_parabola

```
harmonic_basis_from_parabola(
    p: 'Parabola',
    cutoff_energy: 'float'
) -> List[HarmonicOscillator]
```


Create a set of harmonic oscillator wave functions below some cutoff energy.

**p (Parabola)**
The parabola which will be used as the harmonic oscillator's potential surface.

**cutoff_energy (float)**
The energy at which to stop creating basis functions. That is, all basis functions created
will have energy less than or equal to `cutoff_energy`

###get_expansion_coeffs

```
get_expansion_coeffs(state: 'Callable', basis: 'List[Callable]') -> List[float]
```


Given a state and a basis, return the expansion coefficents for that state