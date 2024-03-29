# `ArtraCFD` 

### A Computational Fluid Dynamics Solver
### by 
### Huangrui Mo  ([Google Scholar](https://scholar.google.com/citations?user=zBjRJCgAAAAJ&hl=en)) ([ORCID](https://orcid.org/0000-0001-8279-6329))


## Demo cases

### Case 1: Re=100 flow over a cylinder (click to play animation)

<h1 align="center">
<img width="60%" src="https://github.com/mohuangrui/mohuangrui/blob/main/gallery/1_cyn_vortex_re100.gif" alt="Vortex Shedding">
</h1>

### Case 2: Ma=2.81 Shock diffracting over a cylinder (click to play animation)

<h1 align="center">
<img width="60%" src="https://github.com/mohuangrui/mohuangrui/blob/main/gallery/1_cyn_nomv_invis_m2400_schlieren.gif" alt="Shock Diffraction">
</h1>

### Case 3: Galilean transformation between a supersonic moving wedge and a supersonic incoming flow (click to play animation)

<h1 align="center">
<img width="60%" src="https://github.com/mohuangrui/mohuangrui/blob/main/gallery/1_wedge_deg15_mach2_m1200.gif" alt="Galilean Transformation">
</h1>

### Case 4: Ma=5 flow over a complex body (click to play animation)

<h1 align="center">
<img width="60%" src="https://github.com/mohuangrui/mohuangrui/blob/main/gallery/missile_short_schlieren_m192.gif" alt="Hypersonic Flow">
</h1>

### Case 5: Direct simulation of explosively dispersed particles  (click to play animation)

<h1 align="center">
<img width="60%" src="https://github.com/mohuangrui/mohuangrui/blob/main/gallery/jet_SL2_001_B2_C2_ani.gif" alt="Vortex Shedding">
</h1>

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

## How to compile

1. Download the code
2. Enter the code directory
3. Compile to generate the C executable:
```
make
```

## How to run the program

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

### Reference:

* [Huangrui Mo. Mesoscale modeling and direct simulation of explosively dispersed granular materials. PhD thesis, University of Waterloo, 2019](https://uwspace.uwaterloo.ca/handle/10012/14335)
* [Huangrui Mo, Fue‐Sang Lien, Fan Zhang, and Duane S. Cronin. An immersed boundary method for solving compressible flow with arbitrarily irregular and moving geometry. International Journal for Numerical Methods in Fluids 88, no. 5 (2018): 239-263.](https://onlinelibrary.wiley.com/doi/full/10.1002/fld.4665)
* [Huangrui Mo, Fue-Sang Lien, Fan Zhang, and Duane S. Cronin. A numerical framework for the direct simulation of dense particulate flow under explosive dispersal. Shock Waves 28, no. 3 (2018): 559-577.](https://link.springer.com/article/10.1007/s00193-017-0741-9)
* [Huangrui Mo, Fue-Sang Lien, Fan Zhang, and Duane S. Cronin. A simple field function for solving complex and dynamic fluid-solid system on Cartesian grid. arXiv preprint arXiv:1702.02474 (2017).](https://arxiv.org/abs/1702.02474)
* [Huangrui Mo, Fue-Sang Lien, Fan Zhang, and Duane S. Cronin. A mesoscale study on explosively dispersed granular material using direct simulation. Journal of Applied Physics 125, no. 21 (2019): 214302.](https://aip.scitation.org/doi/full/10.1063/1.5094839)
