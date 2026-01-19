# Ray-Tracing-Code-OCTOPUS
OCTOPUS, a relativistic ray-tracing algorithm developed within a Fortran-based, OpenMP-accelerated framework, designed for asymptotically flat, spherically symmetric curved spacetimes.
This code enables simulations of black hole images, redshift factor images, gravitational lensing effects, light curves, and other related phenomena. 

The code consists of 9 Fortran source files: dynamic.f90, emission.f90, functions.f90, initial_condition.f90, mainprogram.f90, method.f90, model.f90, operations.f90, and parameters.f90.

To implement a new black hole model, users need to:

1. Provide the metric potential and its first three derivatives in the designated section of model.f90

2. Define the corresponding metric parameters in the metric_parameters module within parameters.f90

For complete implementation guidelines and theoretical background, please consult our preprint at arxiv.org (2510.12585):
“OCTOPUS: A Versatile, User-Friendly, and Extensible Public Code for General-Relativistic Ray-Tracing in Spherically Symmetric and Static Spacetimes”.

Release Notes 20260119

1. The previous version employed the Newton method to solve the event horizon equation f(r)=0 with a fixed initial guess of r=2.5. This led to failures (returning NaN) for certain black hole models where a solution could not be found. The current update implements an initial scanning step to bracket the approximate root of f(r)=0. This scan provides a physically reasonable initial guess, after which the Newton iteration proceeds.

2. Fixed an issue where prompt messages or calculation results were not displayed in the terminal window.

3. Added functionality to simulate the trajectory of a single photon, accessible by setting task_model=12. The initial photon coordinates x_start and y_start are defined in the module particle_parameter. The resulting path can be visualized using path_plot.py

4. Resolved an issue that prevented the calculation of the radial velocity for timelike particles within the plunging region of certain spacetime metrics.

