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

 **n (int)** 

> Constraints:  Range(min=0, max=50)

 **center (float)** 

> Constraints:  Range(min=-1000.0, max=1000.0)

 **mass (float)** 

> Default: 1.0

> Constraints:  Range(min=1e-08, max=1000.0)

 **omega (float)** 

> Default: 1.0

> Constraints:  Range(min=1e-08, max=1000.0)



####Properties

#####N


None
#####energy


None
#####potential


None



####Static Methods

#####from\_parabola
```
HarmonicOscillator.from_parabola(
    p: 'Parabola',
    n: 'int',
    mass: 'float' = 1.0
) -> HarmonicOscillator
```


None
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
HarmonicOscillator.__call__(self, x) -> <class 'inspect._empty'>
```


Call self as a function.
#####\_\_kinetic\_\_
```
HarmonicOscillator.__kinetic__(self) -> <class 'inspect._empty'>
```


Return kinetic energy operator applied on this.

K = p^2 / 2m
p = i * sqrt(m w hbar/2)(a+ - a)

k = -1/2 * (m w hbar / 2)[(a+ - a)^2]
[] = (a+^2 + a+a + aa+  - a^2)

[] = sqrt(n+1)sqrt(n+2)| n + 2 >  always
        sqrt(n+1)sqrt(n+1)| n >      always
        sqrt(n)  sqrt(n)  | n >      if n == 0  0
        sqrt(n)  sqrt(n-1)| n - 2 >  if n <= 1  0

k = - (m w hbar) / 4 * []
#####\_\_overlap\_\_
```
HarmonicOscillator.__overlap__(
    self,
    other,
    lower_limit: 'float',
    upper_limit: 'float'
) -> float
```


None

 --- 

###EigenBasis

```
EigenBasis(
    states,
    energies,
    ao_S: 'np.array',
    eigvecs: 'np.array',
    ao_basis: 'List[Callable]'
)
```
#### Fields

 **states (Tuple[Callable, ...])** 

 **energies (Tuple[float, ...])** 

 **ao_S (np.array)** 

 **eigvecs (np.array)** 

 **ao_basis (List[Callable])** 

####Methods

#####transformed
```
EigenBasis.transformed(self, matrix: 'np.array') -> np.array
```


None





####Static Methods

#####from\_basis
```
EigenBasis.from_basis(
    basis: 'List[Callable]',
    H: 'np.array',
    S: 'np.array'
) -> EigenBasis
```


None



####Dunder Methods

#####\_\_len\_\_
```
EigenBasis.__len__(self) -> <class 'inspect._empty'>
```


None



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