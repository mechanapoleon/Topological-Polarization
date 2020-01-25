# Topological-Polarization
Study of the polarization of a linear polymer in a Tight Binding Model, with disordered interaction strenght but fixed position. The study is done for both of open and periodic boundary conditions. In the latter case the polarization can be computed using the Berry phase.

Main code is a Mathematica Notebook, but I have also written some routine in FORTRAN in case I need some faster computation.


## The Physical Problem
We want to compute the polarization of a polymer chain in a tight binding model. The chain is made of an anion and a cation.
For the tight binding model we are consider the on site energy to be + ε and - ε, for the cation and anion respectively. The hopping energy is regulated by a parameter t.

We have three different realization of this physical system
- OBC case: Study the polarization for N cells in open boundary conditions (OBC). We know that for large we should get  P = ±1/2, but we want to know the power law at which it converges.
- Disordered OBC case: same as before, but with the ratio between on site and hopping energy to be taken with a uniform random distribution. We should also look at the energy spectrum: do we get an insulator band energy?
- PBC case: This gets interesting, the polarization is topological and we can compute it by using the Berry phase.


## Theory
We are considering a one dimensional linear lattice, with lattice parameter a. The lattice is made of alternating anions and cations, and we are using a Tight Binding model in order to find the relevant quantities of the system. In brief, the lattice hamiltonian is:



## OBC
In the OBC case we work with localized orbital in position space, the Wannier orbitals. In this program the choice of Wannier orbital is a gaussian. The polarization is easy in this case, it is just straightforwardly the classical operator:

The results we get are:

Also the density of states is interesting:

## Disordered OBC
As we mentioned, for the disordered case we kept the position of the sites fixed but we varied the interaction strenght: the ratio between the on site energy and the hopping energy ε/t is chosen with a uniform random distribution in the interval (0,5). 

Again we have Wannier orbitals, and we compute the same quantities as before.

## PBC
For the periodic boundary conditions we get to the interesting part, which gives name to the repository: the polarization is in fact topological, and we can compute it using the Berry phase.


## References
For the Tight Binding Model:
- Neil W. Ashcroft. N. David Mermin: "Solid State Physics" (1976)
- Giuseppe Grosso, Giuseppe Pastori Parravicini, "Solid State Physics" (2004)

For the application of the Berry phase for computing the polarizat ion
- Nicola A. Spaldin, "A beginner’s guide to the modern theory of polarization", Journal of Solid State Chemistry 195 (2012)
- David Vanderbilt, "Berry Phases in Electronic Structure Theory" (2018)
