# Topological-Polarization
Study of the polarization of a linear polymer in a Tight Binding Model, with disordered interaction strenght but fixed position. The study is done for both of open and periodic boundary conditions. In the latter case the polarization can be computed using the Berry phase.

Main code is a Mathematica Notebook, but I have also written some routine in FORTRAN in case I need some faster computation.


## Problem 
We want to compute the polarization of a polymer chain in a tight binding model. The chain is made of an anion and a cation.
For the tight binding model we are consider the on site energy to be + ε and - ε, for the cation and anion respectively. The hopping energy is regulated by a parameter t.

We have three objectives:
- OBC case: Study the polarization for N cells in open boundary conditions (OBC). We know that for large we should get  P = \pm 1/2, but we are interested how 
- Disordered OBC case: same as before, but with the ratio between on site and hopping energy to be taken with a uniform random distribution. We should also look at the energy spectrum: do we get an insulator band energy?
- PBC case: This gets interesting, the polarization is topological and we can compute it by using the berry phase.


## Theory

## OBC

## Disordered OBC

## PBC

## References
- David Vanderbilt, "Berry Phases in Electronic Structure Theory", Cambridge University Press 2018
