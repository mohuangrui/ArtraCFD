/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_CASEFILEGENERATOR_H_ /* if this is the first definition */
#define ARTRACFD_CASEFILEGENERATOR_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Case Setting File Generator
 *
 * Function
 *      Generate the initial case setting file: artracfd.case and 
 *      initial geometry input file: artracfd.geo for ArtraCFD.
 * Outcome
 *      artracfd.case -- the case setting file.
 *      artracfd.geo -- the initial geometry input file.
 */
extern int GenerateCaseSettingFiles(void);
#endif
/* a good practice: end file with a newline */

