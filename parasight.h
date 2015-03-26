/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_PARASIGHT_H_ /* if this is the first definition */
#define ARTRACFD_PARASIGHT_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Parasight data format and type control, Parasight means the Ensight data
 * format adjusted for Paraview. This adjustment is required for the following
 * reason: the handling of complex geometry in Paraview is achieved by data 
 * filter rather than the blankID of nodes. Therefore, there is no need to
 * write changing geometry files, but need to write a data entry for filter.
 * Because of the same reason, there is no need to write each subparts, since
 * they can be distinguished from each other by the data filter.
 *
 */
typedef char ParasightString[80]; /* Parasight string data requires 80 chars */
typedef float ParasightReal; /* Parasight requires real data to be float */
/*
 * Parasight configuration structure
 */
typedef struct {
    ParasightString baseName; /* data file base name */
    ParasightString fileName; /* store current open file name */
    ParasightString stringData; /* Parasight string data */
}ParasightSet;
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Parasight format data exporter
 *
 * Function
 *      Export conservative field data vector variable:
 *      U = [rho, rho_u, rho_v, rho_w, rho_eT]
 *      to primitive variables = [rho, u, v, w, p, T]
 *      to binary data files with Parasight data format.
 * Notice
 *      U is a linear array that stores all the values.
 *      These data are in sequential state 
 *      and can be accessed by linear index math.
 */
extern int WriteComputedDataParasight(const Real * U, const Space *, 
        const Particle *, const Time *, const Flow *);
/*
 * Parasight format data loader
 *
 * Function
 *      Load computed data from output files which are written in Parasight
 *      format.
 */
extern int LoadComputedDataParasight(Real *U, const Space *, Time *,
        const Flow *);
#endif
/* a good practice: end file with a newline */

