/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_ENSIGHT_H_ /* if this is the first definition */
#define ARTRACFD_ENSIGHT_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Ensight data format and type control
 */
typedef char EnsightString[80]; /* Ensight string data requires 80 chars */
typedef float EnsightReal; /* Ensight requires real data to be float */
/*
 * Ensight configuration structure
 */
typedef struct {
    EnsightString baseName; /* data file base name */
    EnsightString fileName; /* store current open file name */
    EnsightString stringData; /* Ensight string data */
}EnsightSet;
#endif
/* a good practice: end file with a newline */

