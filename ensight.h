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
#include "commons.h"
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
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Ensight transient case file initializer
 *
 * Function
 *      Initialize a Ensight transient case file. This function only needs to
 *      be call once for each non restart run.
 */
extern int InitializeEnsightTransientCaseFile(const Time *);
/*
 * Ensight format data exporter
 *
 * Function
 *      Export conservative field data vector variable:
 *      U = [rho, rho_u, rho_v, rho_w, rho_eT]
 *      to primitive variables = [rho, u, v, w, p, T]
 *      to binary data files with Ensight data format.
 * Notice
 *      U is a linear array that stores all the values.
 *      These data are in sequential state 
 *      and can be accessed by linear index math.
 */
extern int WriteComputedDataEnsight(const Real * U, const Space *, 
        const Particle *, const Time *, const Partition *, const Flow *);
/*
 * Ensight format data loader
 *
 * Function
 *      Load computed data from output files which are written in Ensight
 *      format.
 */
extern int LoadComputedDataEnsight(Real *U, const Space *, Time *,
        const Partition *, const Flow *);
#endif
/* a good practice: end file with a newline */

