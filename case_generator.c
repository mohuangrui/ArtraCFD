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
static int CaseSettingFileGenerator(void);
static int CaseGeometryFileGenerator(void);
static int TriangulatedGeometryFileGenerator(void);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int GenerateCaseFiles(void)
{
    CaseSettingFileGenerator();
    CaseGeometryFileGenerator();
    return 0;
}
static int CaseSettingFileGenerator(void)
{
    FILE *filePointer = fopen("artracfd.case", "w");
    if (NULL == filePointer) {
        FatalError("failed to write case data file: artracfd.case...");
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#                    Case Configuration for ArtraCFD                          -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "# - Coordinate system: Right-handed Cartesian system. X-Y plane is the screen -\n");
    fprintf(filePointer, "#   plane; X is horizontal from west to east; Y is vertical from south to     -\n");
    fprintf(filePointer, "#   north; Z axis is perpendicular to screen and points from front to back.   -\n");
    fprintf(filePointer, "# - Physical quantities are SI Unit based. Data are float type if no specific -\n");
    fprintf(filePointer, "#   information is given. Floats can be exponential notation of 'e'.          -\n");
    fprintf(filePointer, "# - In each 'begin end' environment, there should NOT be any empty or comment -\n");
    fprintf(filePointer, "#   lines. Please double check your input.                                    -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Space Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "space begin\n");
    fprintf(filePointer, "-3, -3, -3         # xmin, ymin, zmin of space domain\n");
    fprintf(filePointer, "3, 3, 3            # xmax, ymax, zmax of space domain (max > min)\n");
    fprintf(filePointer, "250, 250, 1        # x, y, z mesh number (integer; 1: dimension collapse)\n");
    fprintf(filePointer, "space end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Time Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "time begin\n");
    fprintf(filePointer, "0                  # restart number tag (integer; 0: non restart)\n");
    fprintf(filePointer, "1.0                # termination time\n");
    fprintf(filePointer, "1.2                # CFL condition number in (0, 2]\n");
    fprintf(filePointer, "0                  # maximum computing steps (integer; 0: automatic)\n");
    fprintf(filePointer, "1                  # field data writing frequency (integer; 0: infinity)\n");
    fprintf(filePointer, "1                  # data streamer (integer; 0: ParaView; 1: Ensight)\n");
    fprintf(filePointer, "time end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Numerical Method <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# <Type> means the corresponding parameter only takes effect on <Type>\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "numerical begin\n");
    fprintf(filePointer, "1                  # temporal scheme (integer; 0: RK2; 1: RK3;)\n");
    fprintf(filePointer, "1                  # spatial scheme (integer; 0: WENO3; 1: WENO5;)\n");
    fprintf(filePointer, "0                  # multidimensional method (integer; 0: dim split; 1: dim by dim)\n");
    fprintf(filePointer, "0                  # Jacobian average (integer; 0: Arithmetic mean; 1: Roe averages)\n");
    fprintf(filePointer, "0                  # flux splitting method (integer; 0: LLF; 1: SW)\n");
    fprintf(filePointer, "0                  # phase interaction (integer; 0: F; 1: FSI; 2: FSI+SSI)\n");
    fprintf(filePointer, "0                  # layers for reconstruction (integer; 0: infinity)\n");
    fprintf(filePointer, "numerical end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Material Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "material begin\n");
    fprintf(filePointer, "0                  # material (integer; 0: gas; 1: water; 2: solid)\n");
    fprintf(filePointer, "0                  # viscous level (0: none; 1: normal)\n");
    fprintf(filePointer, "0                  # gravity state (integer; 0: off; 1: on)\n");
    fprintf(filePointer, "0, -9.806, 0       # gravity vector\n");
    fprintf(filePointer, "material end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Reference Values  <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "reference begin\n");
    fprintf(filePointer, "1                  # length\n");
    fprintf(filePointer, "1                  # density\n");
    fprintf(filePointer, "1                  # velocity\n");
    fprintf(filePointer, "1                  # temperature\n");
    fprintf(filePointer, "reference end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                             >> Note <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# Values in following parts are relative to reference values. Quantities should\n");
    fprintf(filePointer, "# to be normalized by the given reference values. Like pressure should be\n");
    fprintf(filePointer, "# normalized by reference density times reference velocity square.\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                         >> Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "initialization begin\n");
    fprintf(filePointer, "1                  # density\n");
    fprintf(filePointer, "0                  # x velocity\n");
    fprintf(filePointer, "0                  # y velocity\n");
    fprintf(filePointer, "0                  # z velocity\n");
    fprintf(filePointer, "1                  # pressure\n");
    fprintf(filePointer, "initialization end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Boundary Condition <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# Available types: [inflow], [outflow], [slip wall], [noslip wall], [periodic]\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#west boundary begin\n");
    fprintf(filePointer, "#inflow            # boundary type\n");
    fprintf(filePointer, "#1                 # density\n");
    fprintf(filePointer, "#1                 # x velocity\n");
    fprintf(filePointer, "#0                 # y velocity\n");
    fprintf(filePointer, "#0                 # z velocity\n");
    fprintf(filePointer, "#1                 # pressure\n");
    fprintf(filePointer, "#west boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "west boundary begin\n");
    fprintf(filePointer, "outflow            # boundary type\n");
    fprintf(filePointer, "west boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "east boundary begin\n");
    fprintf(filePointer, "outflow            # boundary type\n");
    fprintf(filePointer, "east boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "south boundary begin\n");
    fprintf(filePointer, "slip wall          # boundary type\n");
    fprintf(filePointer, "-1                 # temperature, negative if adiabatic\n");
    fprintf(filePointer, "south boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "north boundary begin\n");
    fprintf(filePointer, "slip wall          # boundary type\n");
    fprintf(filePointer, "-1                 # temperature, negative if adiabatic\n");
    fprintf(filePointer, "north boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "front boundary begin\n");
    fprintf(filePointer, "outflow            # boundary type\n");
    fprintf(filePointer, "front boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "back boundary begin\n");
    fprintf(filePointer, "outflow            # boundary type\n");
    fprintf(filePointer, "back boundary end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Regional Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#Available options:\n");
    fprintf(filePointer, "#[plane]:          to region on the direction of normal vector\n");
    fprintf(filePointer, "#[sphere]:         to region in the sphere\n");
    fprintf(filePointer, "#[box]:            to region in the box\n");
    fprintf(filePointer, "#[cylinder]:       to region in the cylinder\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "plane initialization begin\n");
    fprintf(filePointer, "-1, 0, 0           # x, y, z of a plane point\n");
    fprintf(filePointer, "-1, 0, 0           # normal vector of plane\n");
    fprintf(filePointer, "3.67372            # density\n");
    fprintf(filePointer, "2.41981            # x velocity\n");
    fprintf(filePointer, "0                  # y velocity\n");
    fprintf(filePointer, "0                  # z velocity\n");
    fprintf(filePointer, "9.04545            # pressure\n");
    fprintf(filePointer, "plane initialization end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "plane initialization begin\n");
    fprintf(filePointer, "-1, 0, 0           # x, y, z of a plane point\n");
    fprintf(filePointer, "1, 0, 0            # normal vector of plane\n");
    fprintf(filePointer, "1                  # density\n");
    fprintf(filePointer, "0                  # x velocity\n");
    fprintf(filePointer, "0                  # y velocity\n");
    fprintf(filePointer, "0                  # z velocity\n");
    fprintf(filePointer, "1                  # pressure\n");
    fprintf(filePointer, "plane initialization end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#sphere initialization begin\n");
    fprintf(filePointer, "#0, 0, 0           # x, y, z of sphere center\n");
    fprintf(filePointer, "#0.1               # radius of sphere\n");
    fprintf(filePointer, "#1                 # density\n");
    fprintf(filePointer, "#0                 # x velocity\n");
    fprintf(filePointer, "#0                 # y velocity\n");
    fprintf(filePointer, "#0                 # z velocity\n");
    fprintf(filePointer, "#1000              # pressure\n");
    fprintf(filePointer, "#sphere initialization end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#box initialization begin\n");
    fprintf(filePointer, "#0, 0, 0           # xmin, ymin, zmin of box\n");
    fprintf(filePointer, "#0.1, 0.1, 0.1     # xmax, ymax, zmax of box\n");
    fprintf(filePointer, "#1                 # density\n");
    fprintf(filePointer, "#0                 # x velocity\n");
    fprintf(filePointer, "#0                 # y velocity\n");
    fprintf(filePointer, "#0                 # z velocity\n");
    fprintf(filePointer, "#1000              # pressure\n");
    fprintf(filePointer, "#box initialization end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#cylinder initialization begin\n");
    fprintf(filePointer, "#0, 0, -0.2        # x1, y1, z1 of center\n");
    fprintf(filePointer, "#0, 0, 0.2         # x2, y2, z2 of center\n");
    fprintf(filePointer, "#0.1               # radius of cylinder\n");
    fprintf(filePointer, "#1                 # density\n");
    fprintf(filePointer, "#0                 # x velocity\n");
    fprintf(filePointer, "#0                 # y velocity\n");
    fprintf(filePointer, "#0                 # z velocity\n");
    fprintf(filePointer, "#1000              # pressure\n");
    fprintf(filePointer, "#cylinder initialization end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Field Data Probes <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "probe count begin\n");
    fprintf(filePointer, "2                  # point probe count (integer; 0: off)\n");
    fprintf(filePointer, "1                  # line probe count (integer; 0: off)\n");
    fprintf(filePointer, "1                  # body-conformal probe (integer; 0: off; 1: on)\n");
    fprintf(filePointer, "1                  # surface force probe (integer; 0: off; 1: on)\n");
    fprintf(filePointer, "probe count end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "probe control begin\n");
    fprintf(filePointer, "1                  # point probe writing frequency (integer; 0: infinity)\n");
    fprintf(filePointer, "1                  # line probe writing frequency (integer; 0: infinity)\n");
    fprintf(filePointer, "1                  # body-conformal probe writing frequency (integer; 0: infinity)\n");
    fprintf(filePointer, "1                  # surface force writing frequency (integer; 0: infinity)\n");
    fprintf(filePointer, "probe control end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "point probe begin\n");
    fprintf(filePointer, "0, -0.5, 0         # x, y, z of the point\n");
    fprintf(filePointer, "0, 0.5, 0          # x, y, z of the point\n");
    fprintf(filePointer, "point probe end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "line probe begin\n");
    fprintf(filePointer, "-0.272, 0.419, 0   # x, y, z of the first end point\n");
    fprintf(filePointer, "2.5, 2.2197, 0     # x, y, z of the second end point\n");
    fprintf(filePointer, "500                # resolution (points on line)\n");
    fprintf(filePointer, "line probe end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#/* a good practice: end file with a newline */\n");
    fprintf(filePointer, "\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int CaseGeometryFileGenerator(void)
{
    FILE *filePointer = fopen("artracfd.geo", "w");
    if (NULL == filePointer) {
        FatalError("failed to write case geometry file: artracfd.geo...");
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#                            Geometry Configuration                           -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "# - Coordinate system: Right-handed Cartesian system. X-Y plane is the screen -\n");
    fprintf(filePointer, "#   plane; X is horizontal from west to east; Y is vertical from south to     -\n");
    fprintf(filePointer, "#   north; Z axis is perpendicular to screen and points from front to back.   -\n");
    fprintf(filePointer, "# - Coordinates and physical quantities should be normalized to dimensionless.-\n");
    fprintf(filePointer, "# - In each 'begin end' environment, there should NOT be any empty or comment -\n");
    fprintf(filePointer, "#   lines. Please double check your input.                                    -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Number of Geometries <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "count begin\n");
    fprintf(filePointer, "1                  # analytical sphere (integer)\n");
    fprintf(filePointer, "1                  # triangulated polyhedron (integer)\n");
    fprintf(filePointer, "count end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                   >> Geometry Information <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# State begin\n");
    fprintf(filePointer, "# O, r, V, W, rho, T, cf, area, volume, mid\n");
    fprintf(filePointer, "# at, ar, ate, g, are, to\n");
    fprintf(filePointer, "# State end\n");
    fprintf(filePointer, "# (Ox, Oy, Oz): geometric center; relative frame for transformation\n");
    fprintf(filePointer, "# r: bounding sphere\n");
    fprintf(filePointer, "# (Vx, Vy, Vz): translational velocity of geometric center\n");
    fprintf(filePointer, "# (Wx, Wy, Wz): rotational velocity relative to geometric center\n");
    fprintf(filePointer, "# rho: density; > 1.0e10 if ignore surface force effect\n");
    fprintf(filePointer, "# T: wall temperature; < 0 if adiabatic; >= 0 if constant\n");
    fprintf(filePointer, "# cf: roughness; <= 0 if slip; > 0 if no slip\n");
    fprintf(filePointer, "# area, volume, mid: surface area, volume, material identifier\n");
    fprintf(filePointer, "# (atx, aty, atz): translational acceleration\n");
    fprintf(filePointer, "# (arx, ary, arz): rotational acceleration\n");
    fprintf(filePointer, "# (atex, atey, atez): exerted external translational acceleration\n");
    fprintf(filePointer, "# (gx, gy, gz): gravitational acceleration\n");
    fprintf(filePointer, "# (arex, arey, arez): exerted external rotational acceleration\n");
    fprintf(filePointer, "# to: time to end external ate and are; <= 0 if never end\n");
    fprintf(filePointer, "# stationary object: V = 0; W = 0; ate = 0; g = 0; are = 0; rho > 1.0e36;\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                 >> Analytical Sphere Section <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "sphere state begin\n");
    fprintf(filePointer, "0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 2700, -1, 1, 0, 0, 0\n");
    fprintf(filePointer, "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n");
    fprintf(filePointer, "sphere state end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                 >> Triangulated Polyhedron Section <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#Polyhedron representation is consistently employed for describing irregular   \n");
    fprintf(filePointer, "#objects. For problems with a collapsed dimension, polyhedron is generated via \n");
    fprintf(filePointer, "#extending the polygon on the collapsed dimension with unit thickness.         \n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "polyhedron geometry begin\n");
    fprintf(filePointer, "artracfd.stl       # geometry file name\n");
    fprintf(filePointer, "polyhedron geometry end\n");
    fprintf(filePointer, "polyhedron state begin\n");
    fprintf(filePointer, "0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 2700, -1, 1, 0, 0, 0\n");
    fprintf(filePointer, "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n");
    fprintf(filePointer, "polyhedron state end\n");
    fprintf(filePointer, "polyhedron transform begin\n");
    fprintf(filePointer, "1, 1, 1, 0, 0, 0, 0, 0, 0 # scale, rotate, translate\n");
    fprintf(filePointer, "polyhedron transform end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#/* a good practice: end file with a newline */\n");
    fprintf(filePointer, "\n");
    fclose(filePointer); /* close current opened file */
    TriangulatedGeometryFileGenerator();
    return 0;
}
static int TriangulatedGeometryFileGenerator(void)
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
    return 0;
}
/* a good practice: end file with a newline */

