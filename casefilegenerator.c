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
#include <string.h> /* manipulating strings */
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
    fprintf(filePointer, "# - Please double check your input!                                           -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "system begin\n");
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
    fprintf(filePointer, "#                        >> Fluid Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "fluid begin\n");
    fprintf(filePointer, "1.205              # density of fluid\n");
    fprintf(filePointer, "1.0e-5             # kinematic viscosity\n");
    fprintf(filePointer, "1.0e2              # thermal diffusivity\n");
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
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Boundary Condition <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "system end\n");
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
    fprintf(filePointer, "# - Please double check your input!                                           -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "geometry begin\n");
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
    fprintf(filePointer, "circle begin\n");
    fprintf(filePointer, "3, 3, 0, 0.5, 0, 0, 0   # x, y, z, radius, u, v, w\n");
    fprintf(filePointer, "circle end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "geometry end\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#/* a good practice: end file with a newline */\n");
    fprintf(filePointer, "\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

