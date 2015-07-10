/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
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
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int CaseSettingFileGenerator(void);
static int CaseGeometryFileGenerator(void);
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
    fprintf(filePointer, "#   information. Floats can be exponential notation of lower-case 'e'.        -\n");
    fprintf(filePointer, "# - In each 'begin end' environment, there should NOT be any empty or comment -\n");
    fprintf(filePointer, "#   lines. Please double check your input!                                    -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Space Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "space begin\n");
    fprintf(filePointer, "-0.5, -0.5, -0.5   # xmin, ymin, zmin of space domain\n");
    fprintf(filePointer, "0.5, 0.5, 0.5      # xmax, ymax, zmax of space domain (max > min)\n");
    fprintf(filePointer, "500, 500, 500      # x, y, z mesh number (integer, 1 if dimension collapse)\n");
    fprintf(filePointer, "space end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Time Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "time begin\n");
    fprintf(filePointer, "0                  # restart flag (integer, 0 if false and 1 if true)\n");
    fprintf(filePointer, "1.0                # total evolution time\n");
    fprintf(filePointer, "-1                 # maximum steps to force cease (integer, -1 if disable)\n");
    fprintf(filePointer, "0.6                # CFL condition number\n");
    fprintf(filePointer, "1                  # total number of times of exporting computed data (integer)\n");
    fprintf(filePointer, "2                  # data streamer (0: ParaView; 1: Ensight; 2 ParaView Ensight)\n");
    fprintf(filePointer, "time end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Numerical Method <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# <Type> means the corresponding parameter only takes effect on <Type>\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "numerical begin\n");
    fprintf(filePointer, "1                  # spatial scheme (0: 2nd Upwind TVD; 1: 5th WENO)\n");
    fprintf(filePointer, "0                  # average method (0: Arithmetic mean; 1: Roe averages)\n");
    fprintf(filePointer, "0                  # <WENO> flux splitting method (0: L-F; 1: S-W)\n");
    fprintf(filePointer, "0.125              # <TVD> Harten's numerical dissipation coefficient [0, 0.5]\n");
    fprintf(filePointer, "numerical end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Fluid and Flow Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "fluid begin\n");
    fprintf(filePointer, "0.71               # Prandtl number\n");
    fprintf(filePointer, "1                  # modify coefficient of dynamic viscosity (0 if inviscid)\n");
    fprintf(filePointer, "fluid end\n");
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
    fprintf(filePointer, "#                             >> NOTE <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# Values in following parts are relative to reference values. Hence, they need\n");
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
    fprintf(filePointer, "periodic           # boundary type\n");
    fprintf(filePointer, "front boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "back boundary begin\n");
    fprintf(filePointer, "periodic           # boundary type\n");
    fprintf(filePointer, "back boundary end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Regional Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#Available options:\n");
    fprintf(filePointer, "#[plane]:          to region on the direction of normal vector points to\n");
    fprintf(filePointer, "#[sphere]:         to region in the sphere\n");
    fprintf(filePointer, "#[box]:            to region in the box\n");
    fprintf(filePointer, "#[cylinder]:       to region in the cylinder\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#plane initialization begin\n");
    fprintf(filePointer, "#0, 0, 0           # x, y, z of a plane point\n");
    fprintf(filePointer, "#-1, 0, 0          # normal vector of plane\n");
    fprintf(filePointer, "#1                 # density\n");
    fprintf(filePointer, "#0                 # x velocity\n");
    fprintf(filePointer, "#0                 # y velocity\n");
    fprintf(filePointer, "#0                 # z velocity\n");
    fprintf(filePointer, "#1000              # pressure\n");
    fprintf(filePointer, "#plane initialization end\n");
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
    fprintf(filePointer, "#Probes originally are all line segments, but can be changed to a single point\n");
    fprintf(filePointer, "#or any allowed number of points by adjusting resolution (points on line).\n");
    fprintf(filePointer, "#NOTICE: the number of probe specify sections should not exceed 10.\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#probe control begin\n");
    fprintf(filePointer, "#1                 # total number of times of exporting probe data (integer)\n");
    fprintf(filePointer, "#probe control end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#probe begin\n");
    fprintf(filePointer, "#-1, 0, 0          # x, y, z of the first end point of line\n");
    fprintf(filePointer, "#1, 0, 0           # x, y, z of the second end point of line\n");
    fprintf(filePointer, "#500               # number of points on line\n");
    fprintf(filePointer, "#probe end\n");
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
    fprintf(filePointer, "#                   Interior Geometry Configuration                           -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "# - Coordinate system: Right-handed Cartesian system. X-Y plane is the screen -\n");
    fprintf(filePointer, "#   plane; X is horizontal from west to east; Y is vertical from south to     -\n");
    fprintf(filePointer, "#   north; Z axis is perpendicular to screen and points from front to back.   -\n");
    fprintf(filePointer, "# - All coordinate values and velocity values here are relative values, they  -\n");
    fprintf(filePointer, "#   should be normalized by the reference length and velocity respectively.   -\n");
    fprintf(filePointer, "# - In each 'begin end' environment, there should NOT be any empty or comment -\n");
    fprintf(filePointer, "#   lines. Please double check your input!                                    -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                   >> Total Number of Objects <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "count begin\n");
    fprintf(filePointer, "1            # total number of objects (integer)\n");
    fprintf(filePointer, "count end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Geometry Information <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "sphere begin\n");
    fprintf(filePointer, "0, 0, 0, 0.5, 1.0e250, 0, 0, 0   # x, y, z, r, rho, u, v, w\n");
    fprintf(filePointer, "sphere end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#/* a good practice: end file with a newline */\n");
    fprintf(filePointer, "\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

