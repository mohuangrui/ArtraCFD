/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_PARAVIEWSTREAM_H_ /* if this is the first definition */
#define ARTRACFD_PARAVIEWSTREAM_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Paraview format data exporter
 *
 * Function
 *      Export conservative field data vector variable:
 *      U = [rho, rho_u, rho_v, rho_w, rho_eT]
 *      to primitive variables = [rho, u, v, w, p, T]
 *      to binary data files with Paraview data format.
 * Notice
 *      U is a linear array that stores all the values.
 *      These data are in sequential state 
 *      and can be accessed by linear index math.
 */
extern int WriteComputedDataParaview(const Real * U, const Space *, 
        const Time *, const Partition *, const Flow *);
/*
 * Paraview data loader
 *
 * Function
 *      Load computed data from output files which are written in Paraview
 *      format.
 */
extern int LoadComputedDataParaview(Real *U, const Space *, Time *,
        const Partition *, const Flow *);
#endif
/* a good practice: end file with a newline */

