/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_PARAVIEW_H_ /* if this is the first definition */
#define ARTRACFD_PARAVIEW_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Paraview data format and type control
 */
typedef char ParaviewString[80]; /* Paraview string data */
typedef double ParaviewReal; /* Paraview real data */
/*
 * Paraview configuration structure
 */
typedef struct {
    ParaviewString baseName; /* data file base name */
    ParaviewString fileName; /* store current open file name */
    ParaviewString floatType; /* Paraview data type */
    ParaviewString byteOrder; /* byte order of data */
}ParaviewSet;
#endif
/* a good practice: end file with a newline */

