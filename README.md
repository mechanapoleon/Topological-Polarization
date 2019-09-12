# Topological-Polarization
Study of the polarization of a linear polymer in a Tight Binding Model, with disordered interaction strenght but fixed position. The study is done for both of open and periodic boundary conditions. In the latter case the polarization can be computed using the Berry phase.

Main code is a Mathematica Notebook, but I have also written some routine in FORTRAN in case I need some faster computation.


## Problem 
We want to compute the polarization of a polymer chain in a tight binding model. The chain is made of an anion and a cation.
For the tight binding model we are consider the on site energy to be + ε and - ε, for the cation and anion respectively. The hopping energy is regulated by a parameter t.

We have three objectives:
- OBC case: Study the polarization for N cells in open boundary conditions (OBC). We know that for large we should get  P = ±1/2, but we want to know the power law at which it converges.
- Disordered OBC case: same as before, but with the ratio between on site and hopping energy to be taken with a uniform random distribution. We should also look at the energy spectrum: do we get an insulator band energy?
- PBC case: This gets interesting, the polarization is topological and we can compute it by using the Berry phase.


## Theory
We are considering a one dimensional linear lattice, with the 

## OBC
In t

## Disordered OBC

## PBC

## References
For the Tight Binding Model:
- Neil W. Ashcroft. N. David Mermin: "Solid State Physics" (1976)
- Giuseppe Grosso, Giuseppe Pastori Parravicini, "Solid State Physics" (2004)

For the application of the Berry phase for computing the polarizat ion
- Nicola A. Spaldin, "A beginner’s guide to the modern theory of polarization", Journal of Solid State Chemistry 195 (2012)
- David Vanderbilt, "Berry Phases in Electronic Structure Theory" (2018)
