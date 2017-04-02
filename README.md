# ArtraCFD
A Computational Fluid Dynamics Solver For Multimaterial Flows

## Solver configuration

### Fluid-solid interaction: Operator splitting

### Fluid dynamics:

#### Governing equations: three-dimensional Navier-Stokes equations (Cartesian, compressible, conservative)
#### Temporal discretization: RK2 and RK3
#### Spatial discretization: WENO3 and WENO5 (convective fluxes) + 2nd order central scheme (diffusive fluxes)
#### Boudary treatment: a sharp interface immersed boundary method (arXiv:1602.06830)

### Solid dynamics:

#### Governing equations: Newton’s second law of translational motion, Euler equations of rotational motion, multi-body contact and collision model
#### Temporal integration: RK2
#### Interface description: triangulated facets with front tracking
