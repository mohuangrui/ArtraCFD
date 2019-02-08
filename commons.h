/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                        *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_COMMONS_H_ /* if undefined */
#define ARTRACFD_COMMONS_H_ /* set a unique marker */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include <stdio.h> /* standard library for input and output */
#include <stdarg.h> /* variable-length argument lists */
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Global integer constants related to governing equation and discretization
 */
typedef enum {
    /* dimensions related to spatial operator */
    DIMO = 4, /* spatial operator dimension: x, y, z operators + source term */
    X = 0,
    Y = 1,
    Z = 2,
    DIMS = 3, /* space dimension: x, y, z */
    PHI = 4, /* source term */
    /* dimension collapse tag */
    COLLAPSEN = 0,
    COLLAPSEX = 1,
    COLLAPSEY = 2,
    COLLAPSEZ = 3,
    COLLAPSEXY = 5,
    COLLAPSEXZ = 7,
    COLLAPSEYZ = 8,
    COLLAPSEXYZ = 17,
    /* dimensions related to temporal operator */
    DIMT = 3, /* number of time levels to store field data */
    TO = 0, /* the time level for current */
    TN = 1, /* the time level for intermediate */
    TM = 2, /* the time level for intermediate */
    /* dimensions related to field variables */
    DIMU = 5, /* conservative vector: rho, rho_u, rho_v, rho_w, rho_eT */
    DIMUo = 6, /* primitive vector: rho, u, v, w, [p, hT, h], [T, c] */
    /* parameters related to numerical model */
    PATHN = 30, /* neighbour searching path */
    PATHSEP = 4, /* layer separator in neighbour searching path: pathN, l1N, l2N, l3N */
    NONE = -1, /* invalid flag */
    WENOTHREE = 0, /* 3rd order weno */
    WENOFIVE = 1, /* 5th order weno */
    OPTSPLIT = 0, /* operator splitting approximation */
    OPTBYOPT = 1, /* operator-by-operator approximation */
    /* parameters related to domain partitions */
    NPART = 15, /* inner region, [west, east, south, north, front, back] x [Boundary, Ghost], physical region, all region */
    PIO = 0, /* the partition region for data iostream */
    PIN = 0,
    PWB = 1,
    PEB = 2,
    PSB = 3,
    PNB = 4,
    PFB = 5,
    PBB = 6,
    PWG = 7,
    PEG = 8,
    PSG = 9,
    PNG = 10,
    PFG = 11,
    PBG = 12,
    PHY = 13,
    PAL = 14,
    LIMIT = 2, /* number of limits */
    MIN = 0,
    MAX = 1,
    /* parameters related to domain boundary conditions */
    NBC = 7, /* Interior, [west, east, south, north, front, back] x [Boundary] */
    INFLOW = 0, /* boundary condition identifier */
    OUTFLOW = 1,
    SLIPWALL = 2,
    NOSLIPWALL = 3,
    PERIODIC = 4,
    VARBC = 6, /* specified primitive variables: rho, u, v, w, p, T */
    /* parameters related to global and regional initialization */
    NIC = 10, /* maximum number of initializer to support */
    ICGLOBAL = 0, /* global initializer */
    ICPLANE = 1, /* plane initializer */
    ICSPHERE = 2, /* sphere initializer */
    ICBOX = 3, /* box initializer */
    ICCYLINDER = 4, /* cylinder initializer */
    POSIC = 7, /* initializer position: x1, y1, z1, [x2, Nx], [y2, Ny], [z2, Nz], r */
    VARIC = 5, /* specified primitive variables: rho, u, v, w, p */
    /* parameters related to geometry */
    DIMTK = 2, /* number of time levels to store kinematic data */
    POLYN = 3, /* polygon facet type */
    EVF = 4, /* edge-vertex-face type */
    /* parameters related to data probes */
    NPROBE = 5, /* point, line, curve, force, space probe */
    PROPT = 0,
    PROLN = 1,
    PROCV = 2,
    PROFC = 3,
    PROSD = 4,
    POSLN = 7, /* x1, y1, z1, x2, y2, z2, resolution */
    /* general parameters */
    STR = 200, /* string length */
    VARSTR =100, /* variable expression length */
} ComConst;
/*
 * Universe data type to improve portability and maintenance
 */
typedef double Real; /* real data */
typedef char String[STR]; /* string data */
typedef int IntVec[DIMS]; /* integer type vector */
typedef Real RealVec[DIMS]; /* real type vector */
/*
 * Global real constants declaration
 */
extern const Real PI;
/*
 * Member structures
 */
typedef struct {
    int did; /* domain identifier */
    int fid; /* closest face identifier */
    int lid; /* interfacial layer identifier */
    int gst; /* ghost layer identifier */
    Real U[DIMT][DIMU]; /* field data at each time level */
} Node; /* field data */

typedef struct {
    IntVec m; /* mesh number of spatial dimensions */
    IntVec n; /* node number of spatial dimensions */
    IntVec ng; /* number of ghost node layers of spatial dimensions */
    int gl; /* number of ghost node layers required by numerical scheme */
    int collapse; /* space collapse flag */
    RealVec d; /* mesh size of spatial dimensions */
    RealVec dd; /* reciprocal of mesh sizes */
    Real tinyL; /* smallest length scale established on grid size */
    int ns[NPART][DIMS][LIMIT]; /* decomposition node range for each partition */
    int np[DIMS][DIMS][LIMIT]; /* computational node range with dimension priority */
    int path[PATHN][DIMS]; /* neighbour searching path */
    int pathSep[PATHSEP]; /* layer separator in neighbour searching path */
    int *restrict typeBC; /* boundary type recorder */
    int (*restrict N)[DIMS]; /* outward surface normal of domain boundary */
    Real (*restrict varBC)[VARBC]; /* field values of each boundary */
    int nIC; /* flow initializer pointer and counter */
    int *restrict typeIC; /* flow initializer type recorder */
    Real (*restrict posIC)[POSIC]; /* position values of each initializer */
    char (*restrict varIC)[VARIC][VARSTR]; /* field expression of each initializer */
    Real domain[DIMS][LIMIT]; /* coordinates define the space domain */
    IntVec proc; /* number of processors of spatial dimensions */
    int procN; /* total number of processors */
} Partition; /* domain discretization and partition */

typedef struct {
    RealVec N; /* normal vector */
    RealVec v0; /* vertex */
    RealVec v1; /* vertex */
    RealVec v2; /* vertex */
} Facet; /* polyhedron facet */

typedef struct {
    int gid; /* geometry identifier */
    IntVec N; /* line of impact */
} Collision; /* collision list */

typedef struct {
    int faceN; /* number of faces. <=0 for analytical polyhedron */
    int edgeN; /* number of edges */
    int vertN; /* number of vertices */
    int state; /* dynamic motion indicator */
    int mid; /* material type */
    Real r; /* bounding sphere radius */
    RealVec O; /* centroid */
    Real I[DIMS][DIMS]; /* inertia matrix */
    Real V[DIMTK][DIMS]; /* translational velocity */
    Real W[DIMTK][DIMS]; /* rotational velocity */
    Real at[DIMTK][DIMS]; /* translational acceleration */
    RealVec g; /* gravitational acceleration */
    Real ar[DIMTK][DIMS]; /* rotational acceleration */
    RealVec Fp; /* pressure force */
    RealVec Fv; /* viscous force */
    RealVec Tt; /* total torque */
    Real to; /* time to end power */
    Real rho; /* density */
    Real T; /* wall temperature */
    Real cf; /* roughness */
    Real area; /* area */
    Real volume; /* volume */
    Real box[DIMS][LIMIT]; /* a bounding box of the polyhedron */
    int (*restrict f)[POLYN]; /* face-vertex list */
    Real (*restrict Nf)[DIMS]; /* face normal */
    int (*restrict e)[EVF]; /* edge-vertex-face list */
    Real (*restrict Ne)[DIMS]; /* edge normal */
    Real (*restrict v)[DIMS]; /* vertex list */
    Real (*restrict Nv)[DIMS]; /* vertex normal */
    Facet *facet; /* facet data */
} Polyhedron; /* polyhedron */

typedef struct {
    int totN; /* total number of geometries */
    int sphN; /* number of analytical polyhedrons */
    int stlN; /* number of triangulated polyhedrons */
    int colN; /* colliding list pointer and count */
    Polyhedron *poly; /* geometry list */
    Collision *col; /* collision list */
} Geometry; /* geometry data */

typedef struct {
    Real eos; /* equation of state */
} Material; /* material property database */
/*
 * Manager structures
 * Memory of normal type members will be automatically allocated from stack.
 * Memory of pointer type members should be dynamically allocated from heap.
 */
typedef struct {
    Node *node; /* field data */
    Geometry geo; /* geometry data */
    Partition part; /* domain discretization and partition data */
} Space;

typedef struct {
    int restart; /* restart tag */
    int stepN; /* total number of steps */
    int stepC; /* step number count */
    int dataN[NPROBE]; /* number for each data probe type */
    int dataW[NPROBE]; /* writing frequency for each data probe type */
    int dataStreamer; /* data streamer */
    int dataC; /* data writing count */
    Real end; /* termination time */
    Real now; /* current time recorder */
    Real numCFL; /* CFL number */
    Real (*restrict pp)[DIMS]; /* point probes */
    Real (*restrict lp)[POSLN]; /* line probes */
} Time;

typedef struct {
    int tScheme; /* temporal discretization scheme */
    int sScheme; /* spatial discretization scheme */
    int sL; /* left offset of stencil index */
    int sR; /* right offset of stencil index */
    int multidim; /* multidimensional space method */
    int jacobMean; /* average method for local Jacobian linearization */
    int fluxSplit; /* flux vector splitting method */
    int psi; /* phase interaction type */
    int ibmLayer; /* number of interfacial layers using flow reconstruction */
    int mid; /* material identifier */
    int gState; /* gravity state */
    int sState; /* source state */
    Real refMa; /* reference Mach number */
    Real refMu; /* reference dynamic viscosity */
    Real gamma; /* heat capacity ratio */
    Real gasR; /* specific gas constant */
    Real cv; /* specific heat capacity at constant volume */
    Real refL; /* characteristic length */
    Real refRho; /* characteristic density */
    Real refV;  /*characteristic velocity */
    Real refT; /* characteristic temperature */
    RealVec g; /* gravity vector */
    Material *mat; /* material database */
} Model;

typedef struct {
    char runMode; /* running mode */
    IntVec proc; /* number of processors per dimension */
} Control;
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Command parser
 *
 * Function
 *      Get rid of end of line and information after #.
 *      Replace tabs with spaces.
 *      Retain only one space between two words.
 *      If no other information exists, produce an empty string.
 */
extern int ParseCommand(char *cmdstr);
/*
 * Format parser
 *
 * Function
 *      Adjust the format string according to the type of Real.
 */
extern char *ParseFormat(char *fmt);
/*
 * Fatal error control
 *
 * Function
 *      Show error and then exit. Once the process exits,
 *      the operating system can free all dynamically allocated
 *      memory associated with the process.
 */
extern void ShowError(const char *fmt, ...);
/*
 * Warning control
 *
 * Function
 *      Show warning.
 */
extern void ShowWarning(const char *fmt, ...);
/*
 * Error stream
 *
 * Function
 *      Accept a variable argument list.
 *      Direct information to standard error.
 */
extern void verror(const char *prefix, const char *fmt, va_list args);
/*
 * Show information
 *
 * Function
 *      Direct information to standard output.
 *      For "Session", it prints a line of asterisks.
 */
extern void ShowInfo(const char *fmt, ...);
/*
 * Assign Storage
 *
 * Function
 *      Return the head address of a linear array of dynamically allocated
 *      memory that is initialized to zero.
 */
extern void *AssignStorage(size_t size);
/*
 * Retrieve storage
 *
 * Function
 *      Free dynamically allocated memory pointed by the pointer.
 */
extern void RetrieveStorage(void *pointer);
/*
 * Auxiliary Functions for File Reading
 *
 * Function
 *      Read in lines until find a matched line.
 *      The file pointer points to the next line of the matched line.
 */
extern void ReadInLine(FILE *fp, const char *line);
/*
 * Auxiliary Functions for File Writing
 *
 * Function
 *      Search down from beginning of file until find a matched line.
 *      The file pointer points to the matched line.
 */
extern void WriteToLine(FILE *fp, const char *line);
/*
 * Standard Stream Functions with Checked Return Values
 */
extern FILE *Fopen(const char *fname, const char *mode);
extern void Fread(void *ptr, size_t size, size_t n, FILE *stream);
extern void Fscanf(FILE *stream, const int n, const char *fmt, ...);
extern void Sscanf(const char *str, const int n, const char *fmt, ...);
/*
 * Read data
 *
 * Function
 *      Read a line from stream and then read n data from the line
 *      and store them according to fmt into the locations pointed
 *      by the elements in the variable argument list.
 *      If n is zero, only read in a line.
 *      If n is negative, no read conversion check.
 */
extern void Sread(FILE *stream, const int n, const char *fmt, ...);
#endif
/* a good practice: end file with a newline */

