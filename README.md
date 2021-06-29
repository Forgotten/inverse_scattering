# Inverse Scattering
Self contained, Light-weight code for solving the inverse scattering problem via optimization based methods.

This repo contains all the necessary code to simulate time-harmonic wave propagation using finite differences. 

In addition, we provide a simple l2 misfit against a reference far-field pattern, and its derivatives using adjoint-state methods. 

For the reconstruction the following methods are implemented

- full-wave form inversion for single frequency using fminunc (quasi-newton) as an optimization routine
- assymptotic least squares

## Organization

- ''src'': contains all the wave propagation and misfit sub routines
- ''examples'': contains the scripts showcasing the inversion algorithms

TODO: 
multi-frequency FWI
least squares inversion
Gauss-Newton optimization

