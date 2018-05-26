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
#ifndef ARTRACFD_COMMONS_H_ /* if this is the first definition */
#define ARTRACFD_COMMONS_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include <stdio.h> /* standard library for input and output */
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Define global integer constants for array bounds, identifiers, etc.
 */
typedef enum {
    /* dimensions related to space */
    DIMS = 3, /* space dimension */
    X = 0,
    Y = 1,
    Z = 2,
    COLLAPSEN = 0, /* dimension collapsed tag */
    COLLAPSEX = 1,
    COLLAPSEY = 2,
    COLLAPSEZ = 3,
    COLLAPSEXY = 5,
    COLLAPSEXZ = 7,
    COLLAPSEYZ = 8, 
    COLLAPSEXYZ = 17, 
    /* dimensions related to field variables */
    DIMU = 5, /* conservative vector: rho, rho_u, rho_v, rho_w, rho_eT */
    DIMUo = 6, /* primitive vector: rho, u, v, w, [p, hT, h], [T, c] */
    DIMT = 3, /* number of time levels to store field data */
    TO = 0, /* the time level for current */
    TN = 1, /* the time level for intermediate */
    TM = 2, /* the time level for intermediate */
    /* parameters related to numerical model */
    PATHN = 30, /* neighbour searching path */
    PATHSEP = 4, /* layer separator in neighbour searching path: pathN, l1N, l2N, l3N */
    NONE = -1, /* invalid flag */
    WENOTHREE = 0, /* 3th order weno */
    WENOFIVE = 1, /* 5th order weno */
    /* parameters related to domain partitions */
    NPART = 13, /* inner region, [west, east, south, north, front, back] x [Boundary, Ghost] */
    NPARTWRITE = 1, /* number of partitions to write data out */
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
    ENTRYBC = 6, /* rho, u, v, w, p, T */
    VARBC = 5, /* rho, u, v, w, p */
    /* parameters related to global and regional initialization */
    NIC = 10, /* maximum number of initializer to support */
    ICGLOBAL = 0, /* global initializer */
    ICPLANE = 1, /* plane initializer */
    ICSPHERE = 2, /* sphere initializer */
    ICBOX = 3, /* box initializer */
    ICCYLINDER = 4, /* cylinder initializer */
    ENTRYIC = 12, /* x1, y1, z1, [x2, Nx], [y2, Ny], [z2, Nz], r, primitive variables */
    VARIC = 5, /* primitive variables: rho, u, v, w, p */
    /* parameters related to geometry */
    DIMTK = 2, /* number of time levels to store kinematic data */
    POLYN = 3, /* polygon facet type */
    EVF = 4, /* edge-vertex-face type */
} Constants;
/*
 * Define some universe data type for portability and maintenance.
 */
typedef double Real; /* real data */
typedef char String[400]; /* string data */
typedef int IntVec[DIMS]; /* integer type vector */
typedef Real RealVec[DIMS]; /* real type vector */
/*
 * Define structures for packing compound data
 */
typedef struct {
    int gid; /* geometry identifier */
    int fid; /* closest face identifier */
    int lid; /* interfacial layer identifier */
    int gst; /* ghost layer identifier */
    Real U[DIMT][DIMU]; /* field data at each time level */
} Node;
/*
 * Domain discretization and partition structure
 */
typedef struct {
    IntVec m; /* mesh number of spatial dimensions */
    IntVec n; /* node number of spatial dimensions */
    int ng; /* number of ghost node layers of global domain */
    int gl; /* number of ghost node layers required for numerical scheme */
    int collapse; /* space collapse flag */
    RealVec d; /* mesh size of spatial dimensions */
    RealVec dd; /* reciprocal of mesh sizes */
    Real tinyL; /* smallest length scale established on grid size */
    int ns[NPART][DIMS][LIMIT]; /* decomposition node range for each partition */
    int np[DIMS][DIMS][LIMIT]; /* computational node range with dimension priority */
    int path[PATHN][DIMS]; /* neighbour searching path */
    int pathSep[PATHSEP]; /* layer separator in neighbour searching path */
    int N[NBC][DIMS]; /* outward surface normal of domain boundary */
    int typeBC[NBC]; /* BC type recorder */
    int countIC; /* flow initializer count */
    int typeIC[NIC]; /* flow initializer type recorder */
    Real valueBC[NBC][ENTRYBC]; /* field values of each boundary */
    Real valueIC[NIC][ENTRYIC]; /* field values of each initializer */
    Real domain[DIMS][LIMIT]; /* coordinates define the space domain */
} Partition;
/*
 * Facet structure
 */
typedef struct {
    RealVec N; /* normal vector */
    RealVec v0; /* vertex */
    RealVec v1; /* vertex */
    RealVec v2; /* vertex */
} Facet;
/*
 * Polyhedron structure
 */
typedef struct {
    int faceN; /* number of faces. 0 for analytical sphere */
    int edgeN; /* number of edges */
    int vertN; /* number of vertices */
    int state; /* dynamic motion indicator */
    int mid; /* material type */
    Real r; /* bounding sphere */
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
} Polyhedron;
/*
 * Collision list
 */
typedef struct {
    int gid; /* geometry identifier */
    IntVec N; /* line of impact */
} Collision;
/*
 * Geometry Entities
 */
typedef struct {
    int totN; /* total number of geometries */
    int sphN; /* number of analytical spheres */
    int stlN; /* number of triangulated polyhedrons */
    int colN; /* colliding list pointer and count */
    Polyhedron *poly; /* geometry list */
    Collision *col; /* collision list */
} Geometry;
/*
 * Material properties
 */
typedef struct {
    Real eos; /* equation of state */
} Material;
/*
 * Space domain parameters
 */
typedef struct {
    Node *node; /* field data */
    Geometry geo; /* geometry in space */
    Partition part; /* domain discretization and partition information */
} Space;
/*
 * Time domain parameters
 */
typedef struct {
    int restart; /* restart tag */
    int stepN; /* total number of steps */
    int stepC; /* step number count */
    int writeN; /* field data writing frequency */
    int writeC; /* field data writing count */
    int dataStreamer; /* types of data streamer */
    int pointWriteN; /* point probe writing frequency */
    int lineWriteN; /* line probe writing frequency */
    int curveWriteN; /* body-conformal probe writing frequency */
    int forceWriteN; /* surface force writing frequency */
    int pointProbeN; /* total number of point probes */
    int lineProbeN; /* total number of line probes */
    int curveProbeN; /* body-conformal probe */
    int forceProbeN; /* surface force probe */
    Real end; /* termination time */
    Real now; /* current time recorder */
    Real numCFL; /* CFL number */
    Real (*restrict pp)[DIMS]; /* point probes */
    Real (*restrict lp)[7]; /* line probes */
} Time;
/*
 * Model properties and physics parameters
 */
typedef struct {
    int tScheme; /* temporal discretization scheme */
    int sScheme; /* spatial discretization scheme */
    int multidim; /* multidimensional space method */
    int jacobMean; /* average method for local Jacobian linearization */
    int fluxSplit; /* flux vector splitting method */
    int fsi; /* material interaction trigger */
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
    Material mat; /* material database */
} Model;
/*
 * Program command line arguments and overall control
 */
typedef struct {
    char runMode; /* mode: [i] interact, [s] serial, [t] threaded, [m] mpi, [g] gpu */
    int procN; /* number of processors */
} Control;
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Command line processor
 *
 * Function
 *      Get rid of end of line, and information after #.
 *      Get rid of before and after tabs, replace tabs with a space.
 *      Get rid of before and after spaces, retain only one space in words.
 *      If no other information exists, the lineCommand turns to a NULL string.
 */
extern int CommandLineProcessor(char *lineCommand);
/*
 * Fatal error control
 *
 * Function
 *      Print information and then exit. Once the process exits, the operating
 *      system is able to free all dynamically allocated memory associated with
 *      the process.
 */
extern void FatalError(const char *statement);
/*
 * Show information to terminal
 *
 * Function
 *      Print information to standard out. Statement is the information to show. 
 *      If statement is "Session End", it prints a line asterisks.
 */
extern int ShowInformation(const char *statement);
/*
 * Assign Storage
 *
 * Function
 *      Use malloc to assign a linear array of storage. Returns the head address
 *      of the assigned storage. Since malloc does not initialize the storage, 
 *      a call of memset is used to initialize the assigned memory to zero.
 */
extern void *AssignStorage(size_t size);
/*
 * Retrieve storage
 *
 * Function
 *      Use free to free the storage space of the pointer.
 * Notice
 *      Don't free pointer of storage that not allocated by dynamic allocation.
 *      The original pointer becomes to be a wild pointer after being freed, be
 *      aware of this situation. It's a better practice to set pointer back to 
 *      NULL after calling free.
 */
extern int RetrieveStorage(void *pointer);
/*
 * Auxiliary Functions for File Reading
 *
 * Function
 *      Read in lines from the current line until a line matches the lineString.
 *      The file pointer points to the next line of the matched line.
 */
extern int ReadInLine(FILE *filePointer, const char *lineString);
/*
 * Auxiliary Functions for File Writing
 *
 * Function
 *      Search down the file from beginning until a line matches the lineString.
 *      The file pointer points to the matched line.
 */
extern int WriteToLine(FILE *filePointer, const char *lineString);
/*
 * Standard Stream Functions with Checked Return Values
 */
extern void Fgets(char *str, int num, FILE *stream);
extern void Fread(void *ptr, size_t size, size_t count, FILE *stream);
extern void VerifyReadConversion(const int num, const int expect);
#endif
/* a good practice: end file with a newline */

