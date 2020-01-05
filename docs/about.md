# About

This library is a reproduction of a project I did as a student
of computational chemistry. In essence, the project's goal 
was to find an answer to 

> How long does it take for an electron to cross a molecular bridge?

Students of Quantum Mechanics will soon realize that this is question
it not even really a valid one. In the realm of QM, there is no
such thing as a particle being at "location A". Rather, a particle
will have some probability distribution of existing at any location, 
which is only verifiable upon measurement. When drafting the 
manuscript for this work, my advisor and myself had considerable
difficulty in simply describing the problem in language that wasn't 
flat-out wrong. 


In retrospect, this application is really solving the **hitting time**
of a time-heterogenous markov chain, wherein the transition matrix
is the probability that a discretized chunk of probability density
is the "object" that is hopping between spatial regions. 


This codebase is not the same code that was used in the manuscript, 
but rather a reimagining that I've done several years later, 
after much more experience in software development. I did this 
to see how different the code looks like after some (5 years) of 
experience in creating software, with 2 of those being professionally. 


Working the field of Machine Learning has really opened my eyes
to the power of well engineered libraries in the pursuit of science. 
The deep learning library keras in particular is a triumph of API
design greatly accelerating the progress of research in the field. 


The problems faced in machine learning and those faced in computational
physics and chemistry have an extremely large degree of overlap. While 
generating this library I've tried to create something that would
be intuitive for a domain expert to use. 

```python
    from transit_chem.utils import Parabola
    from transit_chem.basis import harmonic_basis_from_parabola, HarmonicOscillator, EigenBasis
    from transit_chem import potentials
    from transit_chem.time_evolution import TimeEvolvingState
    from transit_chem.operators import Hamiltonian, Overlap, Kinetic


    parabola = Parabola(a=1.0, b=1.0, c=1.0)
    basis = harmonic_basis_from_parabola(parabola, cutoff_energy=4.0)
    v = potentials.Harmonic(center=0.0, mass=1.0, omega=1.0)
    H = Hamiltonian(v).matrix(basis)
    S = Overlap().matrix(basis)
    eig_basis = EigenBasis.from_basis(basis, H, S)
    initial_state = HarmonicOscillator(n=0, center=parabola.vertex)
    time_evolving = TimeEvolvingState(eig_basis, initial_state)
    kinetic_over_time = time_evolving.observable(Kinetic(), hermitian=True)
    print(kinetic_over_time(t=0))
    print(kinetic_over_time(t=10))

```

The above is valid, executable code.
If I've done my job correctly, this code snippet will be easy for
domain experts to understand.