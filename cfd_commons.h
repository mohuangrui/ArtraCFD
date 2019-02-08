/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                        *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_CFD_COMMONS_H_ /* if undefined */
#define ARTRACFD_CFD_COMMONS_H_ /* set a unique marker */
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
 * Average method
 *
 * Function
 *      Compute averaged variables at interface.
 */
extern void SymmetricAverage(const int averager, const Real gamma,
        const Real UL[restrict], const Real UR[restrict], Real Uo[restrict]);
/*
 * Jacobian matrices, eigenvalues, and eigenvectors
 *
 * Function
 *      Compute eigenvalues and eigenvectors.
 */
extern void Eigenvalue(const int s, const Real Uo[restrict], Real Lambda[restrict]);
extern void EigenvalueSplitting(const int splitter, const Real Lambda[restrict],
        Real LambdaP[restrict], Real LambdaN[restrict]);
extern void EigenvectorL(const int s, const Real gamma, const Real Uo[restrict],
        Real L[restrict][DIMU]);
extern void EigenvectorR(const int s, const Real Uo[restrict], Real R[restrict][DIMU]);
/*
 * Convective fluxes
 *
 * Function
 *      Compute convective fluxes.
 */
extern void ConvectiveFlux(const int s, const Real gamma, const Real U[restrict], Real F[restrict]);
/*
 * Physical property
 */
extern Real Viscosity(const Real T);
extern Real PrandtlNumber(void);
/*
 * Compute the values of primitive variable vector
 *
 * Function
 *      Compute primitive variable vector according to conservative vector.
 */
extern void MapPrimitive(const Real gamma, const Real gasR, const Real U[restrict], Real Uo[restrict]);
extern Real ComputePressure(const Real gamma, const Real U[restrict]);
extern Real ComputeTemperature(const Real cv, const Real U[restrict]);
/*
 * Compute and update conservative variable vector
 *
 * Function
 *      Compute conservative variable vector according to primitive vector.
 */
extern void MapConservative(const Real gamma, const Real Uo[restrict], Real U[restrict]);
/*
 * Index math
 *
 * Function
 *      Calculate the node index.
 */
extern int IndexNode(const int k, const int j, const int i, const int jMax, const int iMax);
/*
 * Verify node region
 *
 * Function
 *     Check whether a node is within the part box.
 */
extern int InPartBox(const int k, const int j, const int i, const int pbox[restrict][LIMIT]);
/*
 * Coordinates transformation
 *
 * Function
 *      Transform coordinates between node coordinates and general coordinates.
 */
extern int MapNode(const Real s, const Real sMin, const Real dds, const int ng);
extern int ConfineSpace(const int n, const int nMin, const int nMax);
extern Real MapPoint(const int n, const Real sMin, const Real ds, const int ng);
/*
 * Common math functions
 */
extern Real MinReal(const Real x, const Real y);
extern Real MaxReal(const Real x, const Real y);
extern int EqualReal(const Real x, const Real y);
extern int MinInt(const int x, const int y);
extern int MaxInt(const int x, const int y);
extern int Sign(const Real x);
extern Real Dot(const Real V1[restrict], const Real V2[restrict]);
extern Real Norm(const Real V[restrict]);
extern Real Dist2(const Real V1[restrict], const Real V2[restrict]);
extern Real Dist(const Real V1[restrict], const Real V2[restrict]);
extern void Cross(const Real V1[restrict], const Real V2[restrict], Real V[restrict]);
extern void OrthogonalSpace(const Real N[restrict], Real Ta[restrict], Real Tb[restrict]);
extern void Normalize(const int dimV, const Real normalizer, Real V[restrict]);
#endif
/* a good practice: end file with a newline */

