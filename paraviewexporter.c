/****************************************************************************
 * Export Computed Data in VTK format                                       *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "paraview.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int InitializeParaviewDataFile(ParaviewSet *paraSet, const Time *);
static int WriteParaviewDataFile(ParaviewSet *paraSet, const Time *);
static int WriteParticleFile(ParaviewSet *paraSet, const Particle *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteComputedDataParaview(const Real *U, const Space *space, 
        const Particle *particle, const Time *time, const Partition *part, 
        const Flow *flow)
{
    ShowInformation("  writing field data to file...");
    ParaviewSet paraSet = { /* initialize ParaviewSet environment */
        .baseName = "paraview", /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .stringData = {'\0'}, /* string data recorder */
    };
    if (0 == time->stepCount) { /* this is the initialization step */
        InitializeParaviewDataFile(&paraSet, time);
    }
    WriteParaviewDataFile(&paraSet, time);
    WriteParaviewVariableFile(U, &paraSet, space, flow);
    WriteParticleFile(&paraSet, particle);
    return 0;
}
static int InitializeParaviewDataFile(ParaviewSet *paraSet, const Time *time)
{
    FILE *filePointer = fopen("paraview.pvd", "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to transient case file...");
    }
    /* output information to file */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"Collection\" version=\"0.1\"\n");
    fprintf(filePointer, "         byte_order=\"LittleEndian\"\n");
    fprintf(filePointer, "         compressor=\"vtkZLibDataCompressor\">\n");
    fprintf(filePointer, "  <Collection>\n");
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteParaviewDataFile(ParaviewSet *paraSet, const Time *time)
{
    /* store updated basename in filename */
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%05d", paraSet->baseName, time->outputCount); 
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
    fprintf(filePointer, "<VTKFile type=\"Collection\" version=\"0.1\"\n");
    fprintf(filePointer, "         byte_order=\"LittleEndian\"\n");
    fprintf(filePointer, "         compressor=\"vtkZLibDataCompressor\">\n");
    fprintf(filePointer, "  <Collection>\n");
    fprintf(filePointer, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", time->currentTime);
    fprintf(filePointer, "             file=\"%s.vtp\"/>\n", paraSet->baseName);
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "  <!-- stepCount %d -->\n", time->stepCount);
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
    for (int line = 1; line < targetLine; ++line) {
        fgets(currentLine, sizeof currentLine, filePointer);
    }
    /* append informatiom */
    fprintf(filePointer, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", time->currentTime);
    fprintf(filePointer, "             file=\"%s.vtp\"/>\n", paraSet->baseName);
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteParticleFile(ParaviewSet *paraSet, const Particle *particle)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.particle", paraSet->baseName);
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("faild to write particle data file...");
    }
    fprintf(filePointer, "N: %d\n", particle->totalN); /* number of objects */
    const Real *ptk = NULL;
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = particle->headAddress + geoCount * particle->entryN; /* point to storage of current particle */
        fprintf(filePointer, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                ptk[0], ptk[1], ptk[2], ptk[3], ptk[4], ptk[5], 
                ptk[6], ptk[7]);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

