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
#ifndef ARTRACFD_CFD_COMMONS_H_ /* if this is the first definition */
#define ARTRACFD_CFD_COMMONS_H_ /* a unique marker for this header file */
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
 * Jacobian matrices, eigenvalues, and eigenvectors
 *
 * Function
 *      Compute eigenvalues and eigenvectors.
 */
extern void Eigenvalue(const int s, const Real Uo[], Real Lambda[]);
extern void EigenvalueSplitting(const int splitter, const Real Lambda[], Real LambdaP[], Real LambdaN[]);
extern void EigenvectorL(const int s, const Real gamma, const Real Uo[], Real L[][DIMU]);
extern void EigenvectorR(const int s, const Real Uo[], Real R[][DIMU]);
/*
 * Average method
 *
 * Function
 *      Compute averaged variables at interface.
 */
extern void SymmetricAverage(const int averager, const Real gamma, const Real UL[], const Real UR[], Real Uo[]);
/*
 * Convective fluxes
 *
 * Function
 *      Compute convective fluxes.
 */
extern void ConvectiveFlux(const int s, const Real gamma, const Real U[], Real F[]);
/*
 * Diffusive fluxes
 *
 * Function
 *      Compute Diffusive fluxes.
 */
extern void DiffusiveFlux(const int s, const int k, const int j, const int i, const Space *, const Model *, Real Fv[]);
/*
 * Compute the values of primitive variable vector
 *
 * Parameter
 *      Uo[] -- a array stores the returned values of primitives.
 * Notice
 *      calculated values are [rho, u, v, w, p]
 */
extern void PrimitiveByConservative(const Real gamma, const Real gasR, const Real U[], Real Uo[]);
extern Real ComputePressure(const Real gamma, const Real U[]);
extern Real ComputeTemperature(const Real cv, const Real U[]);
/*
 * Compute and update conservative variable vector
 *
 * Function
 *      Compute and update conservative variable vector according to primitive values.
 */
extern void ConservativeByPrimitive(const Real gamma, const Real Uo[], Real U[]);
/*
 * Index math
 *
 * Function
 *      calculate the index of current node.
 * Returns
 *      int -- the calculated index value
 */
extern int IndexNode(const int k, const int j, const int i, const int jMax, const int iMax);
/*
 * Coordinates transformation
 *
 * Function
 *      transform coordinates between node coordinates and general coordinates.
 */
extern int NodeSpace(const Real sMin, const Real s, const Real dds, const int ng);
extern int ValidNodeSpace(const int ns, const int nsMin, const int nsMax);
extern Real PointSpace(const Real sMin, const int ns, const Real ds, const int ng);
/*
 * Common math functions
 */
extern Real MinReal(const Real x, const Real y);
extern Real MaxReal(const Real x, const Real y);
extern int MinInt(const int x, const int y);
extern int MaxInt(const int x, const int y);
extern int Sign(const Real x);
extern Real Dot(const RealVector V1, const RealVector V2);
extern Real Norm(const RealVector V);
extern Real Dist2(const RealVector V1, const RealVector V2);
extern Real Dist(const RealVector V1, const RealVector V2);
extern void Cross(const RealVector V1, const RealVector V2, RealVector V);
extern void OrthogonalSpace(const RealVector N, RealVector Ta, RealVector Tb);
extern void Normalize(const int dimV, const Real normalizer, Real V[]);
#endif
/* a good practice: end file with a newline */

