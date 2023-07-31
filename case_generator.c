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
 * Required Header Files
 ****************************************************************************/
#include "case_generator.h"
#include <stdio.h> /* standard library for input and output */
#include "stl.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void GenerateCaseSettingFile(void);
static void GenerateCaseGeometryFile(void);
static void GenerateTriangulatedGeometryFile(void);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void GenerateCaseFiles(void)
{
    GenerateCaseSettingFile();
    GenerateCaseGeometryFile();
    return;
}
static void GenerateCaseSettingFile(void)
{
    FILE *fp = Fopen("artracfd.case", "w");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "#                    Case Configuration for ArtraCFD                          -\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "# - Coordinate system: Right-handed Cartesian system. X-Y plane is the screen -\n");
    fprintf(fp, "#   plane; X is horizontal from west to east; Y is vertical from south to     -\n");
    fprintf(fp, "#   north; Z axis is perpendicular to screen and points from front to back.   -\n");
    fprintf(fp, "# - Physical quantities are SI Unit based. Data are float type if no specific -\n");
    fprintf(fp, "#   information is given. Floats can be exponential notation of 'e'.          -\n");
    fprintf(fp, "# - In each 'begin end' environment, there should NOT be any empty or comment -\n");
    fprintf(fp, "#   lines. Please double check input.                                         -\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                          >> Space Domain <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "space begin\n");
    fprintf(fp, "-3, -3, -3         # xmin, ymin, zmin\n");
    fprintf(fp, "3, 3, 3            # xmax, ymax, zmax (max > min)\n");
    fprintf(fp, "250, 250, 1        # mx, my, mz (int; 1: dimension collapse)\n");
    fprintf(fp, "space end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                          >> Time Domain <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "time begin\n");
    fprintf(fp, "0                  # restart data checkpoint (int; 0: none)\n");
    fprintf(fp, "1.0                # termination time\n");
    fprintf(fp, "1.2                # CFL condition number in (0, 2]\n");
    fprintf(fp, "0                  # maximum computing steps (int; 0: auto)\n");
    fprintf(fp, "1                  # space data writing frequency (int; 0: inf)\n");
    fprintf(fp, "1                  # data streamer (int; 0: ParaView; 1: Ensight)\n");
    fprintf(fp, "time end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                        >> Numerical Method <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "numerical begin\n");
    fprintf(fp, "1                  # temporal scheme (int; 0: RK2; 1: RK3;)\n");
    fprintf(fp, "1                  # spatial scheme (int; 0: WENO3; 1: WENO5;)\n");
    fprintf(fp, "0                  # dimension scheme (int; 0: dim split; 1: dim by dim)\n");
    fprintf(fp, "0                  # Jacobian average (int; 0: Arithmetic; 1: Roe)\n");
    fprintf(fp, "0                  # flux splitting method (int; 0: LLF; 1: SW)\n");
    fprintf(fp, "0                  # phase interaction (int; 0: F; 1: FSI; 2: FSI+SSI)\n");
    fprintf(fp, "1                  # ibm reconstruction layers (int; 0: inf)\n");
    fprintf(fp, "numerical end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                        >> Material Properties <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "material begin\n");
    fprintf(fp, "0                  # material type (int; 0: gas; 1: water; 2: solid)\n");
    fprintf(fp, "0                  # viscous level (0: none; 1: normal)\n");
    fprintf(fp, "0                  # gravity state (int; 0: off; 1: on)\n");
    fprintf(fp, "0, -9.806, 0       # gravity vector\n");
    fprintf(fp, "material end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                        >> Reference Values  <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "reference begin\n");
    fprintf(fp, "1                  # length\n");
    fprintf(fp, "1                  # density\n");
    fprintf(fp, "1                  # velocity\n");
    fprintf(fp, "1                  # temperature\n");
    fprintf(fp, "reference end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                             >> Note <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# Physical quantities below should be normalized by the reference values.\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                         >> Initialization <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "initialization begin\n");
    fprintf(fp, "1                  # density expression\n");
    fprintf(fp, "0                  # x velocity expression\n");
    fprintf(fp, "0                  # y velocity expression\n");
    fprintf(fp, "0                  # z velocity expression\n");
    fprintf(fp, "1                  # pressure expression\n");
    fprintf(fp, "initialization end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                        >> Boundary Condition <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# Available types: [inflow], [outflow], [slip wall], [noslip wall], [periodic]\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#west boundary begin\n");
    fprintf(fp, "#inflow            # boundary type\n");
    fprintf(fp, "#1                 # density\n");
    fprintf(fp, "#1                 # x velocity\n");
    fprintf(fp, "#0                 # y velocity\n");
    fprintf(fp, "#0                 # z velocity\n");
    fprintf(fp, "#1                 # pressure\n");
    fprintf(fp, "#west boundary end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "west boundary begin\n");
    fprintf(fp, "outflow            # boundary type\n");
    fprintf(fp, "west boundary end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "east boundary begin\n");
    fprintf(fp, "outflow            # boundary type\n");
    fprintf(fp, "east boundary end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "south boundary begin\n");
    fprintf(fp, "slip wall          # boundary type\n");
    fprintf(fp, "-1                 # temperature (<0: adiabatic)\n");
    fprintf(fp, "south boundary end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "north boundary begin\n");
    fprintf(fp, "slip wall          # boundary type\n");
    fprintf(fp, "-1                 # temperature (<0: adiabatic)\n");
    fprintf(fp, "north boundary end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "front boundary begin\n");
    fprintf(fp, "outflow            # boundary type\n");
    fprintf(fp, "front boundary end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "back boundary begin\n");
    fprintf(fp, "outflow            # boundary type\n");
    fprintf(fp, "back boundary end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                  >> Regional Initialization <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#Available options:\n");
    fprintf(fp, "#[plane]:          to region on the direction of normal vector\n");
    fprintf(fp, "#[sphere]:         to region in the sphere\n");
    fprintf(fp, "#[box]:            to region in the box\n");
    fprintf(fp, "#[cylinder]:       to region in the cylinder\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "plane initialization begin\n");
    fprintf(fp, "-1, 0, 0           # x, y, z of a plane point\n");
    fprintf(fp, "-1, 0, 0           # normal vector of plane\n");
    fprintf(fp, "3.67372            # density expression\n");
    fprintf(fp, "2.41981            # x velocity expression\n");
    fprintf(fp, "0                  # y velocity expression\n");
    fprintf(fp, "0                  # z velocity expression\n");
    fprintf(fp, "9.04545            # pressure expression\n");
    fprintf(fp, "plane initialization end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "plane initialization begin\n");
    fprintf(fp, "-1, 0, 0           # x, y, z of a plane point\n");
    fprintf(fp, "1, 0, 0            # normal vector of plane\n");
    fprintf(fp, "1                  # density expression\n");
    fprintf(fp, "0                  # x velocity expression\n");
    fprintf(fp, "0                  # y velocity expression\n");
    fprintf(fp, "0                  # z velocity expression\n");
    fprintf(fp, "1                  # pressure expression\n");
    fprintf(fp, "plane initialization end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#sphere initialization begin\n");
    fprintf(fp, "#0, 0, 0           # x, y, z of sphere center\n");
    fprintf(fp, "#0.1               # radius of sphere\n");
    fprintf(fp, "#1                 # density expression\n");
    fprintf(fp, "#0                 # x velocity expression\n");
    fprintf(fp, "#0                 # y velocity expression\n");
    fprintf(fp, "#0                 # z velocity expression\n");
    fprintf(fp, "#1000              # pressure expression\n");
    fprintf(fp, "#sphere initialization end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#box initialization begin\n");
    fprintf(fp, "#0, 0, 0           # xmin, ymin, zmin of box\n");
    fprintf(fp, "#0.1, 0.1, 0.1     # xmax, ymax, zmax of box\n");
    fprintf(fp, "#1                 # density expression\n");
    fprintf(fp, "#0                 # x velocity expression\n");
    fprintf(fp, "#0                 # y velocity expression\n");
    fprintf(fp, "#0                 # z velocity expression\n");
    fprintf(fp, "#1000              # pressure expression\n");
    fprintf(fp, "#box initialization end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#cylinder initialization begin\n");
    fprintf(fp, "#0, 0, -0.2        # x1, y1, z1 of center\n");
    fprintf(fp, "#0, 0, 0.2         # x2, y2, z2 of center\n");
    fprintf(fp, "#0.1               # radius of cylinder\n");
    fprintf(fp, "#1                 # density expression\n");
    fprintf(fp, "#0                 # x velocity expression\n");
    fprintf(fp, "#0                 # y velocity expression\n");
    fprintf(fp, "#0                 # z velocity expression\n");
    fprintf(fp, "#1000              # pressure expression\n");
    fprintf(fp, "#cylinder initialization end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                    >> Field Data Probes <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "probe count begin\n");
    fprintf(fp, "2                  # point probe count (int; 0: off)\n");
    fprintf(fp, "1                  # line probe count (int; 0: off)\n");
    fprintf(fp, "1                  # body-conformal probe (int; 0: off; 1: on)\n");
    fprintf(fp, "1                  # surface force probe (int; 0: off; 1: on)\n");
    fprintf(fp, "probe count end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "probe control begin\n");
    fprintf(fp, "1                  # point probe writing frequency (int; 0: inf)\n");
    fprintf(fp, "1                  # line probe writing frequency (int; 0: inf)\n");
    fprintf(fp, "1                  # body-conformal probe writing frequency (int; 0: inf)\n");
    fprintf(fp, "1                  # surface force writing frequency (int; 0: inf)\n");
    fprintf(fp, "probe control end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "point probe begin\n");
    fprintf(fp, "0, -0.5, 0         # x, y, z\n");
    fprintf(fp, "0, 0.5, 0          # x, y, z\n");
    fprintf(fp, "point probe end\n");
    fprintf(fp, "#\n");
    fprintf(fp, "line probe begin\n");
    fprintf(fp, "-0.272, 0.419, 0   # x1, y1, z1\n");
    fprintf(fp, "2.5, 2.2197, 0     # x2, y2, z2\n");
    fprintf(fp, "500                # resolution\n");
    fprintf(fp, "line probe end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#/* a good practice: end file with a newline */\n");
    fprintf(fp, "\n");
    fclose(fp);
    return;
}
static void GenerateCaseGeometryFile(void)
{
    FILE *fp = Fopen("artracfd.geo", "w");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "#                            Geometry Configuration                           -\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "# - Coordinate system: Right-handed Cartesian system. X-Y plane is the screen -\n");
    fprintf(fp, "#   plane; X is horizontal from west to east; Y is vertical from south to     -\n");
    fprintf(fp, "#   north; Z axis is perpendicular to screen and points from front to back.   -\n");
    fprintf(fp, "# - Coordinates and physical quantities should be normalized to dimensionless.-\n");
    fprintf(fp, "# - In each 'begin end' environment, there should NOT be any empty or comment -\n");
    fprintf(fp, "#   lines. Please double check input.                                         -\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                  >> Number of Geometries <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "count begin\n");
    fprintf(fp, "1                  # analytical polyhedron (int)\n");
    fprintf(fp, "1                  # triangulated polyhedron (int)\n");
    fprintf(fp, "count end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                   >> Geometry Information <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "# State begin\n");
    fprintf(fp, "# O, r, V, W, rho, T, cf, area, volume, mid\n");
    fprintf(fp, "# at, ar, ate, g, are, to\n");
    fprintf(fp, "# State end\n");
    fprintf(fp, "# (Ox, Oy, Oz): geometric center; relative frame for transformation\n");
    fprintf(fp, "# r: bounding sphere radius\n");
    fprintf(fp, "# (Vx, Vy, Vz): translational velocity of geometric center\n");
    fprintf(fp, "# (Wx, Wy, Wz): rotational velocity relative to geometric center\n");
    fprintf(fp, "# rho: density; > 1.0e10 if ignore surface force effect\n");
    fprintf(fp, "# T: wall temperature; < 0 if adiabatic; >= 0 if constant\n");
    fprintf(fp, "# cf: roughness; <= 0 if slip; > 0 if no slip\n");
    fprintf(fp, "# area, volume, mid: surface area, volume, material identifier\n");
    fprintf(fp, "# (atx, aty, atz): translational acceleration\n");
    fprintf(fp, "# (arx, ary, arz): rotational acceleration\n");
    fprintf(fp, "# (atex, atey, atez): exerted external translational acceleration\n");
    fprintf(fp, "# (gx, gy, gz): gravitational acceleration\n");
    fprintf(fp, "# (arex, arey, arez): exerted external rotational acceleration\n");
    fprintf(fp, "# to: time to end external ate and are; <= 0 if never end\n");
    fprintf(fp, "# stationary object: V = 0; W = 0; ate = 0; g = 0; are = 0; rho > 1.0e36;\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                 >> Analytical Polyhedron Section <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "sphere state begin\n");
    fprintf(fp, "0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 2700, -1, 1, 0, 0, 0\n");
    fprintf(fp, "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n");
    fprintf(fp, "sphere state end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                 >> Triangulated Polyhedron Section <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#Polyhedron representation is consistently employed for describing irregular   \n");
    fprintf(fp, "#objects. For problems with a collapsed dimension, polyhedron is generated via \n");
    fprintf(fp, "#extending the polygon on the collapsed dimension with unit thickness.         \n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "polyhedron geometry begin\n");
    fprintf(fp, "artracfd.stl       # geometry file name\n");
    fprintf(fp, "polyhedron geometry end\n");
    fprintf(fp, "polyhedron state begin\n");
    fprintf(fp, "0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 2700, -1, 1, 0, 0, 0\n");
    fprintf(fp, "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n");
    fprintf(fp, "polyhedron state end\n");
    fprintf(fp, "polyhedron transform begin\n");
    fprintf(fp, "1, 1, 1, 0, 0, 0, 0, 0, 0 # scale, rotate, translate\n");
    fprintf(fp, "polyhedron transform end\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#/* a good practice: end file with a newline */\n");
    fprintf(fp, "\n");
    fclose(fp);
    GenerateTriangulatedGeometryFile();
    return;
}
static void GenerateTriangulatedGeometryFile(void)
{
    Facet facetData[8] = {
        {{-5.000000e-001, 8.660254e-001, 0.000000e+000},
            {8.660254e-001, 5.000000e-001, 5.000000e-001},
            {8.660254e-001, 5.000000e-001, -5.000000e-001},
            {0.000000e+000, 0.000000e+000, 5.000000e-001}},
        {{-5.000000e-001, 8.660254e-001, 0.000000e+000},
            {0.000000e+000, 0.000000e+000, 5.000000e-001},
            {8.660254e-001, 5.000000e-001, -5.000000e-001},
            {0.000000e+000, 0.000000e+000, -5.000000e-001}},
        {{1.000000e+000, 0.000000e+000, 0.000000e+000},
            {8.660254e-001, -5.000000e-001, 5.000000e-001},
            {8.660254e-001, -5.000000e-001, -5.000000e-001},
            {8.660254e-001, 5.000000e-001, 5.000000e-001}},
        {{1.000000e+000, 0.000000e+000, 0.000000e+000},
            {8.660254e-001, 5.000000e-001, 5.000000e-001},
            {8.660254e-001, -5.000000e-001, -5.000000e-001},
            {8.660254e-001, 5.000000e-001, -5.000000e-001}},
        {{-5.000000e-001, -8.660254e-001, 0.000000e+000},
            {0.000000e+000, 0.000000e+000, 5.000000e-001},
            {0.000000e+000, 0.000000e+000, -5.000000e-001},
            {8.660254e-001, -5.000000e-001, 5.000000e-001}},
        {{-5.000000e-001, -8.660254e-001, 0.000000e+000},
            {8.660254e-001, -5.000000e-001, 5.000000e-001},
            {0.000000e+000, 0.000000e+000, -5.000000e-001},
            {8.660254e-001, -5.000000e-001, -5.000000e-001}},
        {{0.000000e+000, 0.000000e+000, -1.000000e+000},
            {8.660254e-001, 5.000000e-001, -5.000000e-001},
            {8.660254e-001, -5.000000e-001, -5.000000e-001},
            {0.000000e+000, 0.000000e+000, -5.000000e-001}},
        {{0.000000e+000, 0.000000e+000, 1.000000e+000},
            {0.000000e+000, 0.000000e+000, 5.000000e-001},
            {8.660254e-001, -5.000000e-001, 5.000000e-001},
            {8.660254e-001, 5.000000e-001, 5.000000e-001}}
    };
    Polyhedron wedge = {
        .faceN = 8,
        .facet = facetData
    };
    WriteStlFile("artracfd.stl", &wedge);
    return;
}
/* a good practice: end file with a newline */

