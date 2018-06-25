# `ArtraCFD` A Computational Fluid Dynamics Solver

## Solver configuration

### Fluid-solid interaction:

* Operator splitting

### Fluid dynamics:

* Governing equations: 3D Navier-Stokes equations (Cartesian, compressible, conservative)
* Temporal discretization: RK2 and RK3
* Spatial discretization: WENO3 and WENO5 (convective fluxes) + 2nd order central scheme (diffusive fluxes)
* Boudary treatment: immersed boundary method

### Solid dynamics:

* Governing equations: Newton's second law (translation), Euler equations (rotation), multi-body contact and collision
* Temporal integration: RK2
* Interface description: triangulated facets with front tracking

