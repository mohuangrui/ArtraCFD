/****************************************************************************
 * Geometry Data Loader and Writer                                          *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a loader and writer for geometry data                *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "geometrystream.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "paraview.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NonrestartGeometryLoader(Particle *);
static int RestartGeometryLoader(Particle *);
static int LoadGeometryDataParaview(Particle *);
static int AllocateMemoryForGeometryData(Particle *);
static int ComputeGeometryParameters(Particle *, const Space *, const Flow *);
static int WriteGeometryDataParaview(const Particle *, const Time *);
static int InitializeTransientParaviewDataFile(ParaviewSet *);
static int WriteSteadyParaviewDataFile(ParaviewSet *, const Time *);
static int WriteParaviewVariableFile(const Particle *, ParaviewSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function load geometry data from the geometry file.
 */
int LoadGeometryData(Particle *particle, const Space *space, const Time *time,
        const Flow *flow)
{
    if (0 == time->restart) { /* if non-restart, read input file */
        NonrestartGeometryLoader(particle);
    } else { /* if restart, read the particle file */
        RestartGeometryLoader(particle);
    }
    ComputeGeometryParameters(particle, space, flow);
    return 0;
}
static int NonrestartGeometryLoader(Particle *particle)
{
    ShowInformation("Loading geometry data ...");
    FILE *filePointer = fopen("artracfd.geo", "r");
    if (NULL == filePointer) {
        FatalError("failed to open geometry file: artracfd.geo...");
    }
    /* read and process file line by line */
    char currentLine[500] = {'\0'}; /* store the current read line */
    Real *ptk = NULL;
    int entryCount = 0; /* entry count */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "sphere count begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(particle->totalN)); 
            continue;
        }
        if (0 == strncmp(currentLine, "sphere begin", sizeof currentLine)) {
            ++entryCount;
            if (0 == particle->totalN) { /* no internal geometries */
                continue;
            }
            AllocateMemoryForGeometryData(particle);
            /* set format specifier according to the type of Real */
            char format[100] = "%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg"; /* default is double type */
            if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
                strncpy(format, "%g, %g, %g, %g, %g, %g, %g, %g", sizeof format); /* float type */
            }
            for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
                ptk = IndexGeometry(geoCount, particle) ;
                fgets(currentLine, sizeof currentLine, filePointer);
                sscanf(currentLine, format, ptk + 0, ptk + 1, ptk + 2, ptk + 3, ptk + 4, ptk + 5,
                        ptk + 6, ptk + 7);
            }
        }
        continue;
    }
    fclose(filePointer); /* close current opened file */
    /* Check missing information section in configuration */
    if (2 != entryCount) {
        FatalError("missing or repeated necessary information section");
    }
    ShowInformation("Session End");
    return 0;
}
static int AllocateMemoryForGeometryData(Particle *particle)
{
    /* 
     * Assign storage to store particle information:
     * x, y, z, r, density, u, v, w,       fx, fy, fz, tally  area  1/mass
     * 0, 1, 2, 3,    4,    5, 6, 7,        8,  9, 10,  11     12     13           ENTRYPTK: 14
     *    need to be read in                  calculated
     */
    particle->headAddress = AssignStorage(particle->totalN * ENTRYPTK, "Real");
    return 0;
}
static int RestartGeometryLoader(Particle *particle)
{
    ShowInformation("Restore geometry data ...");
    LoadGeometryDataParaview(particle);
    ShowInformation("Session End");
    return 0;
}
static int LoadGeometryDataParaview(Particle *particle)
{
    FILE *filePointer = NULL;
    Real *ptk = NULL;
    filePointer = fopen("restart.geo", "r");
    if (NULL == filePointer) {
        FatalError("failed to reload geometry from restart.geo file...");
    }
    /* get rid of redundant lines */
    char currentLine[200] = {'\0'}; /* store current line */
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    /* get total number of particles */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %d", &(particle->totalN)); 
    if (0 == particle->totalN) { /* no internal geometries */
        fclose(filePointer); /* close current opened file */
        return 0;
    }
    AllocateMemoryForGeometryData(particle);
    ParaviewReal data = 0.0; /* the Paraview data format */
    /* set format specifier according to the type of Real */
    char format[5] = "%lg"; /* default is double type */
    if (sizeof(ParaviewReal) == sizeof(float)) {
        strncpy(format, "%g", sizeof format); /* float type */
    }
    for (int dim = 0; dim < 8; ++dim) {
        fgets(currentLine, sizeof currentLine, filePointer);
        for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
            fscanf(filePointer, format, &data);
            ptk = IndexGeometry(geoCount, particle);
            ptk[dim] = data;
        }
        fgets(currentLine, sizeof currentLine, filePointer); /* get rid of the end of line of data */
        fgets(currentLine, sizeof currentLine, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int ComputeGeometryParameters(Particle *particle, const Space *space, const Flow *flow)
{
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        Real *ptk = IndexGeometry(geoCount, particle);
        if (1 == space->collapsed) { /* space dimension collapsed */
            ptk[12] = 2.0 * ptk[3] * flow->pi; /* circle perimeter */
            ptk[13] = ptk[4] * ptk[3] * ptk[3] * flow->pi; /* circle mass */
        } else {
            ptk[12] = 4.0 * ptk[3] * ptk[3] * flow->pi; /* sphere surface */
            ptk[13] = ptk[4] * (4.0 / 3.0) * ptk[3] * ptk[3] * ptk[3] * flow->pi; /* sphere mass */
        }
        /* get the mass reciprocal */
        ptk[13] = 1 / ptk[13];
    }
    return 0;
}
/*
 * Write geometry information of particles.
 */
int WriteGeometryData(const Particle *particle, const Time *time)
{
    WriteGeometryDataParaview(particle, time);
    return 0;
}
static int WriteGeometryDataParaview(const Particle *particle, const Time *time)
{
    ParaviewSet paraSet = { /* initialize ParaviewSet environment */
        .baseName = "particle", /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .floatType = "Float32", /* paraview data type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    if (0 == time->stepCount) { /* this is the initialization step */
        InitializeTransientParaviewDataFile(&paraSet);
    }
    WriteSteadyParaviewDataFile(&paraSet, time);
    WriteParaviewVariableFile(particle, &paraSet);
    return 0;
}
static int InitializeTransientParaviewDataFile(ParaviewSet *paraSet)
{
    FILE *filePointer = fopen("particle.pvd", "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to transient case file...");
    }
    /* output information to file */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"Collection\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <Collection>\n");
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteSteadyParaviewDataFile(ParaviewSet *paraSet, const Time *time)
{
    /* store updated basename in filename */
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%05d", 
            paraSet->baseName, time->outputCount); 
    /* basename is updated here! */
    snprintf(paraSet->baseName, sizeof(ParaviewString), "%s", paraSet->fileName); 
    /* current filename */
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.pvd", paraSet->baseName); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to steady case file...");
    }
    /* output information to file */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"Collection\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <Collection>\n");
    fprintf(filePointer, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", 
            time->currentTime);
    fprintf(filePointer, "             file=\"%s.vts\"/>\n", paraSet->baseName);
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "  <!-- Order %d -->\n", time->outputCount);
    fprintf(filePointer, "  <!-- Time %.6g -->\n", time->currentTime);
    fprintf(filePointer, "  <!-- Step %d -->\n", time->stepCount);
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    /*
     * Add the current export to the transient case
     */
    filePointer = fopen("paraview.pvd", "r+");
    if (NULL == filePointer) {
        FatalError("failed to add data to transient file...");
    }
    /* seek the target line for adding information */
    char currentLine[200] = {'\0'}; /* store the current read line */
    int targetLine = 1;
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "</Collection>", sizeof currentLine)) {
            break;
        }
        ++targetLine;
    }
    /* redirect to the target line */
    rewind(filePointer); /* seek to the beginning of the file */
    for (int line = 1; line < targetLine; ++line) {
        fgets(currentLine, sizeof currentLine, filePointer);
    }
    /* append informatiom */
    fprintf(filePointer, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", 
            time->currentTime);
    fprintf(filePointer, "             file=\"%s.vts\"/>\n", paraSet->baseName);
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteParaviewVariableFile(const Particle *particle, 
        ParaviewSet *paraSet)
{
    FILE *filePointer = NULL;
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.vts", paraSet->baseName); 
    filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to write data file...");
    }
    Real *ptk = NULL;
    ParaviewReal data = 0.0; /* paraview scalar data */
    ParaviewReal vector[3] = {0.0}; /* paraview vector data elements */
    /* the scalar values at each node in current part */
    const char name[12][5] = {"x", "y", "z", "r", "rho", "u", "v", "w", "fx", "fy", "fz", "id"};
    int iMin = 0;
    int iMax = particle->totalN - 1;
    int jMin = 0;
    int jMax = 0;
    int kMin = 0;
    int kMax = 0;
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"StructuredGrid\" version=\"0.1\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 
            iMin, iMax, jMin, jMax, kMin, kMax);
    fprintf(filePointer, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", 
            iMin, iMax, jMin, jMax, kMin, kMax);
    fprintf(filePointer, "      <PointData Scalars=\"r\" Vectors=\"Vel\">\n");
    fprintf(filePointer, "      <!-- N %d -->\n", particle->totalN);
    for (int dim = 0; dim < 12; ++dim) {
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"%s\" format=\"ascii\">\n", 
                paraSet->floatType, name[dim]);
        fprintf(filePointer, "          ");
        for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
            ptk = IndexGeometry(geoCount, particle);
            if (11 == dim) {
                data = geoCount;
            } else {
                data = ptk[dim];
            }
            fprintf(filePointer, "%.6g ", data);
        }
        fprintf(filePointer, "\n        </DataArray>\n");
    }
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"Vel\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = IndexGeometry(geoCount, particle);
        vector[0] = ptk[5];
        vector[1] = ptk[6];
        vector[2] = ptk[7];
        fprintf(filePointer, "%.6g %.6g %.6g ", vector[0], vector[1], vector[2]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"Force\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = IndexGeometry(geoCount, particle);
        vector[0] = ptk[8];
        vector[1] = ptk[9];
        vector[2] = ptk[10];
        fprintf(filePointer, "%.6g %.6g %.6g ", vector[0], vector[1], vector[2]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </PointData>\n");
    fprintf(filePointer, "      <CellData>\n");
    fprintf(filePointer, "      </CellData>\n");
    fprintf(filePointer, "      <Points>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"points\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = IndexGeometry(geoCount, particle);
        vector[0] = ptk[0];
        vector[1] = ptk[1];
        vector[2] = ptk[2];
        fprintf(filePointer, "%.6g %.6g %.6g ", vector[0], vector[1], vector[2]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </Points>\n");
    fprintf(filePointer, "    </Piece>\n");
    fprintf(filePointer, "  </StructuredGrid>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

