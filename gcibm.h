/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_GCIBM_H_ /* if this is the first definition */
#define ARTRACFD_GCIBM_H_ /* a unique marker for this header file */
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
 * Compute domain geometry
 *
 * Function
 *      Employ GCIBM approach to handle complex geometry that locates in
 *      the computational domain.
 */
extern int ComputeDomainGeometryGCIBM(Space *, Particle *, const Partition *);
/*
 * Boundary condition for interior ghost cells
 *
 * Function
 *      Apply boundary conditions to the interior ghost cells.
 */
extern int BoundaryConditionGCIBM(Real *U, const Space *, const Particle *, 
        const Partition *, const Flow *);
/*
 * Inverse Distance Weighting
 *
 * Function
 *      Reconstruction of the values of primitive vector Uo for a spatial
 *      point (z, y, x) based on the neighbours around node (k, j, i) in
 *      index range h by inversed distance weighting.
 *
 * Return
 *      Weighted values without normalization, the normalization factor is
 *      saved and returned at the last element of Uo.
 */
extern int InverseDistanceWeighting(Real Uo[], const Real z, const Real y, const Real x,
        const int k, const int j, const int i, const int h, const Real *U, 
        const Space *, const Flow *);
/*
 * In geometry criteria
 *
 * Function
 *      Check node whether locates in the geometry pointed by the geometry
 *      pointer.
 *
 * Return
 *      negative -- in current geometry
 *      zero     -- on current geometry
 *      positive -- out of current geometry
 */
extern Real InGeometry(const int k, const int j, const int i, const Real *ptk, const Space *);
/*
 * Calculate geometry information
 *
 * Function
 *
 *      Calculate geometry information of current node.
 */
extern int CalculateGeometryInformation(Real info[], const int k, const int j, 
        const int i, const Real *ptk, const Space *);
#endif
/* a good practice: end file with a newline */

