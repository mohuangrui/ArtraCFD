# `ArtraCFD` A Computational Fluid Dynamics Solver

## Solver configuration

### Fluid-solid interaction:

* Operator splitting

### Fluid dynamics:

* Governing equations: three-dimensional Navier-Stokes equations (Cartesian, compressible, conservative)
* Temporal discretization: RK2 and RK3
* Spatial discretization: WENO3 and WENO5 (convective fluxes) + 2nd order central scheme (diffusive fluxes)
* Boudary treatment: immersed boundary method

### Solid dynamics:

* Governing equations: Newton’s second law (translation), Euler equations (rotation), multi-body contact and collision
* Temporal integration: RK2
* Interface description: triangulated facets with front tracking

### References:

* Mo, H., Lien, F. S., Zhang, F., & Cronin, D. S. (2017). A numerical framework for the direct simulation of dense particulate flow under explosive dispersal. Shock Waves, 28(3) 559-577.
* Mo, H., Lien, F. S., Zhang, F., & Cronin, D. S. (2018). A simple field function for solving complex and dynamic fluid-solid system on Cartesian grid. arXiv preprint arXiv:1702.02474.
* Mo, H., Lien, F. S., Zhang, F., & Cronin, D. S. (2018). An immersed boundary method for solving compressible flow with arbitrarily irregular and moving geometry. arXiv preprint arXiv:1602.06830.
 
