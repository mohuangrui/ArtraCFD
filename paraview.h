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
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Paraview data format and type control
 */
typedef char ParaviewString[80]; /* Paraview string data */
typedef float ParaviewReal; /* Paraview real data */
/*
 * Paraview configuration structure
 */
typedef struct {
    ParaviewString baseName; /* data file base name */
    ParaviewString fileName; /* store current open file name */
    ParaviewString stringData; /* Paraview string data */
}ParaviewSet;
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
        const Particle *, const Time *, const Flow *);
/*
 * Paraview data loader
 *
 * Function
 *      Load computed data from output files which are written in VTK
 *      format.
 */
extern int LoadComputedDataParaview(Real *U, const Space *, Time *,
        const Flow *);
#endif
/* a good practice: end file with a newline */

