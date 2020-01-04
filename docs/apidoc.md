# transit_chem

I haven't figured these API docs out. You're better off reading
the source code from the github repo. 


# transit_chem.basis

## HarmonicOscillator
```python
HarmonicOscillator(self, n: 'int', center: 'float', mass: 'float' = 1.0, omega: 'float' = 1.0) -> None
```
A 1D quantum harmonic oscillator wave function.

Parameters
-----------
n
    The quantum number
center
    The center of the potential well
omega
    The angular frequency of the oscillator
mass
    The mass of the particle

Examples
---------

>>> ho = HarmonicOscillator(n=2, center=0.5)
>>> ho
HarmonicOscillator(n=2, center=0.5, omega=1, mass=1.0)
>>> ho(0.5)
-0.5311259660135984
>>> ho(1000)
0.0
>>> ho(-1000)
0.0


## get_expansion_coeffs
```python
get_expansion_coeffs(state: 'Callable', basis: 'List[Callable]') -> 'List[float]'
```
Memoized version of get_expansion_coeffs(state: 'Callable', basis: 'List[Callable]') -> 'List[float]'


