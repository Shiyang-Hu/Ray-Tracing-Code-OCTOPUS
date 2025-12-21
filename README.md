# Ray-Tracing-Code-OCTOPUS
OCTOPUS, a relativistic ray-tracing algorithm developed within a Fortran-based, OpenMP-accelerated framework, designed for asymptotically flat, spherically symmetric curved spacetimes.
This code enables simulations of black hole images, redshift factor images, gravitational lensing effects, light curves, and other related phenomena. 

The code consists of 9 Fortran source files: dynamic.f90, emission.f90, functions.f90, initial_condition.f90, mainprogram.f90, method.f90, model.f90, operations.f90, and parameters.f90.

To implement a new black hole model, users need to:

1. Provide the metric potential and its first three derivatives in the designated section of model.f90

2. Define the corresponding metric parameters in the metric_parameters module within parameters.f90

For complete implementation guidelines and theoretical background, please consult our preprint at arxiv.org (2510.12585):
“OCTOPUS: A Versatile, User-Friendly, and Extensible Public Code for General-Relativistic Ray-Tracing in Spherically Symmetric and Static Spacetimes”.
