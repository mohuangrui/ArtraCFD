/****************************************************************************
 * Generate and Initialize Case Setting Files of ArtraCFD                   *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - Generates the input files of artracfd: artracfd.case and artracfd.geo  *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "casefilegenerator.h"
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
int GenerateCaseSettingFiles(void)
{
    CaseSettingFileGenerator();
    CaseGeometryFileGenerator();
    return 0;
}
static int CaseSettingFileGenerator(void)
{
    FILE *filePointer = fopen("artracfd.case", "w");
    if (filePointer == NULL) {
        FatalError("failed to write case data file: artracfd.case...");
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#                    Case Configuration for ArtraCFD                          -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "# - Coordinate system: Right-handed Cartesian system. X-Y plane is the screen -\n");
    fprintf(filePointer, "#   plane; X is horizontal from west to east; Y is vertical from south to     -\n");
    fprintf(filePointer, "#   north; Z axis is perpendicular to the screen and points from front to     -\n");
    fprintf(filePointer, "#   back; The origin locates at the west-south-front corner of the            -\n");
    fprintf(filePointer, "#   computational domain;                                                     -\n");
    fprintf(filePointer, "# - Physical quantities are SI Unit based. Data are float type if no specific -\n");
    fprintf(filePointer, "#   information. Floats can be exponential notation of lower-case e.          -\n");
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
    fprintf(filePointer, "6, 6, 0            # x, y, z length (0 if collapse a dimension)\n");
    fprintf(filePointer, "100, 100, 1        # x, y, z mesh number (integer, 1 if collapse a dimension)\n");
    fprintf(filePointer, "2                  # number of exterior ghost cells layers (integer, >=1)\n");
    fprintf(filePointer, "space end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Time Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "time begin\n");
    fprintf(filePointer, "0                  # restart flag (integer, 0 if false and 1 if true)\n");
    fprintf(filePointer, "1.4                # total evolution time\n");
    fprintf(filePointer, "0.4                # CFL condition number\n");
    fprintf(filePointer, "2                  # total number of times of exporting computed data (integer)\n");
    fprintf(filePointer, "time end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Fluid and Flow Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "fluid begin\n");
    fprintf(filePointer, "0.71               # Prandtl number\n");
    fprintf(filePointer, "fluid end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Reference Values  <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "reference begin\n");
    fprintf(filePointer, "1                  # length\n");
    fprintf(filePointer, "1.205              # density\n");
    fprintf(filePointer, "340                # velocity\n");
    fprintf(filePointer, "293                # temperature\n");
    fprintf(filePointer, "reference end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                            >> NOTE <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# Values in following parts are relative to reference values. Hence, they need\n");
    fprintf(filePointer, "# to be normalized by the given reference values. Like pressure should be\n");
    fprintf(filePointer, "# normalized by reference density times reference velocity square.\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Flow Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "initialization begin\n");
    fprintf(filePointer, "1                 # density\n");
    fprintf(filePointer, "1                 # x velocity\n");
    fprintf(filePointer, "0                 # y velocity\n");
    fprintf(filePointer, "0                 # z velocity\n");
    fprintf(filePointer, "1                 # pressure\n");
    fprintf(filePointer, "initialization end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Boundary Condition <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# Available boundary types: [inlet], [outflow], [slip wall], [nonslip wall]\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "west boundary begin\n");
    fprintf(filePointer, "inlet             # boundary type\n");
    fprintf(filePointer, "1                 # density\n");
    fprintf(filePointer, "1                 # x velocity\n");
    fprintf(filePointer, "0                 # y velocity\n");
    fprintf(filePointer, "0                 # z velocity\n");
    fprintf(filePointer, "1                 # pressure\n");
    fprintf(filePointer, "west boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "east boundary begin\n");
    fprintf(filePointer, "outflow           # boundary type\n");
    fprintf(filePointer, "east boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "south boundary begin\n");
    fprintf(filePointer, "slip wall         # boundary type\n");
    fprintf(filePointer, "-1                # temperature, negative if adiabatic\n");
    fprintf(filePointer, "south boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "north boundary begin\n");
    fprintf(filePointer, "slip wall         # boundary type\n");
    fprintf(filePointer, "-1                # temperature, negative if adiabatic\n");
    fprintf(filePointer, "north boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "front boundary begin\n");
    fprintf(filePointer, "slip wall         # boundary type\n");
    fprintf(filePointer, "-1                # temperature, negative if adiabatic\n");
    fprintf(filePointer, "front boundary end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "back boundary begin\n");
    fprintf(filePointer, "slip wall         # boundary type\n");
    fprintf(filePointer, "-1                # temperature, negative if adiabatic\n");
    fprintf(filePointer, "back boundary end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Special Regional Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#Available options:\n");
    fprintf(filePointer, "#[plane]:  to region on the plane or the direction of normal vector points to\n");
    fprintf(filePointer, "#[sphere]: to region on the sphere or in the sphere\n");
    fprintf(filePointer, "#[box]:    to region on the box or in the box\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#plane initialization begin\n");
    fprintf(filePointer, "#3, 3, 3           # x, y, z of a plane point\n");
    fprintf(filePointer, "#1, 1, 1           # normal vector of plane\n");
    fprintf(filePointer, "#1                 # density\n");
    fprintf(filePointer, "#1                 # x velocity\n");
    fprintf(filePointer, "#0                 # y velocity\n");
    fprintf(filePointer, "#0                 # z velocity\n");
    fprintf(filePointer, "#1                 # pressure\n");
    fprintf(filePointer, "#plane initialization end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#sphere initialization begin\n");
    fprintf(filePointer, "#3, 3, 3           # x, y, z of sphere center\n");
    fprintf(filePointer, "#1                 # radius of sphere\n");
    fprintf(filePointer, "#1                 # density\n");
    fprintf(filePointer, "#1                 # x velocity\n");
    fprintf(filePointer, "#0                 # y velocity\n");
    fprintf(filePointer, "#0                 # z velocity\n");
    fprintf(filePointer, "#1                 # pressure\n");
    fprintf(filePointer, "#sphere initialization end\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#box initialization begin\n");
    fprintf(filePointer, "#2, 2, 2           # x, y, z of west-south-front corner of box\n");
    fprintf(filePointer, "#4, 4, 4           # x, y, z of east-north-back corner of box\n");
    fprintf(filePointer, "#1                 # density\n");
    fprintf(filePointer, "#1                 # x velocity\n");
    fprintf(filePointer, "#0                 # y velocity\n");
    fprintf(filePointer, "#0                 # z velocity\n");
    fprintf(filePointer, "#1                 # pressure\n");
    fprintf(filePointer, "#box initialization end\n");
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
    if (filePointer == NULL) {
        FatalError("failed to write case geometry file: artracfd.geo...");
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#                   Internal Geometry Configuration                           -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "# - Coordinate system: Right-handed Cartesian system. X-Y plane is the screen -\n");
    fprintf(filePointer, "#   plane; X is horizontal from west to east; Y is vertical from south to     -\n");
    fprintf(filePointer, "#   north; Z axis is perpendicular to the screen and points from front to     -\n");
    fprintf(filePointer, "#   back; The origin locates at the west-south-front corner of the            -\n");
    fprintf(filePointer, "#   computational domain;                                                     -\n");
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
    fprintf(filePointer, "sphere count begin\n");
    fprintf(filePointer, "1            # total number of objects (integer)\n");
    fprintf(filePointer, "sphere count end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Geometry Information <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "shpere begin\n");
    fprintf(filePointer, "3, 3, 0, 0.5, 0, 0, 0   # x, y, z, radius, u, v, w\n");
    fprintf(filePointer, "sphere end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#/* a good practice: end file with a newline */\n");
    fprintf(filePointer, "\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

