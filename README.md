# `ArtraCFD` A Computational Fluid Dynamics Solver

### By [Huangrui Mo](https://orcid.org/0000-0001-8279-6329)

### How to compile

1. Download the code
2. Enter the code directory
3. Compile to generate the C executable:
```
make
```

### How to run the program

1. Run the program:
```
./artracfd
```
2. Check help guide: 
```
help
```
3. Check program manual:
```
manual
```
4. Generate a sample test case:
```
init
```
5. Solve the smaple case:
```
solve
```

For algorithms and more test cases, please check the `Reference` below.

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

### Reference:

* [Huangrui Mo. Mesoscale modeling and direct simulation of explosively dispersed granular materials. PhD thesis, University of Waterloo, 2019](https://uwspace.uwaterloo.ca/handle/10012/14335)
