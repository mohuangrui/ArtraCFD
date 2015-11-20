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
extern int DecompositionCoefficientAlpha(Real alpha[], Real L[][DIMU],
        const int idx, const int idxh, const Real *U);
extern int EigenvalueLambda(const int s, Real lambda[], const Real Uo[]);
extern int EigenvectorSpaceL(const int s, Real L[][DIMU], const Real Uo[], const Real gamma);
extern int EigenvectorSpaceR(const int s, Real R[][DIMU], const Real Uo[]);
/*
 * Flux splitting
 *
 * Function
 *      Splitting eigenvalues by flux splitting methods.
 */
extern int FluxSplitting(Real lambdaPlus[], Real lambdaMinus[], const Real lambda[], const int splitter);
/*
 * Average method
 *
 * Function
 *      Compute averaged variables at interface.
 */
extern int SymmetricAverage(Real Uo[], const int idx, const int idxh,
        const Real *U, const Real gamma, const int averager);
/*
 * Convective fluxes
 *
 * Function
 *      Compute convective fluxes.
 */
extern int ConvectiveFlux(const int s, Real F[], const int idx, const Real *U, const Real gamma);
/*
 * Diffusive fluxes
 *
 * Function
 *      Compute Diffusive fluxes.
 */
extern int DiffusiveFlux(const int s, Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *, const Model *);
/*
 * Compute the values of primitive variable vector
 *
 * Parameter
 *      Uo[] -- a array stores the returned values of primitives.
 *      idx  -- the index of current node.
 * Notice
 *      calculated values are [rho, u, v, w, p]
 */
extern int PrimitiveByConservative(Real Uo[], const int idx, const Real *U, const Model *);
extern Real ComputePressure(const int idx, const Real *U, const Model *);
extern Real ComputeTemperature(const int idx, const Real *U, const Model *);
/*
 * Compute and update conservative variable vector
 *
 * Function
 *      Compute and update conservative variable vector according to primitive values.
 */
extern int ConservativeByPrimitive(Real *U, const int idx, const Real Uo[], const Model *);
/*
 * Index math
 *
 * Function
 *      calculate the index of current node.
 * Returns
 *      int -- the calculated index value
 */
extern int IndexNode(const int k, const int j, const int i, const Space *);
/*
 * Coordinates transformation
 *
 * Function
 *      transform coordinates between node coordinates and general coordinates.
 */
extern int NodeSpace(const Real point, const int s, const Space *);
extern int ValidNodeSpace(const int node, const int s, const Partition *);
extern Real PointSpace(const int node, const int s, const Space *);
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
extern int Cross(RealVector V, const RealVector V1, const RealVector V2);
extern int OrthogonalSpace(const RealVector N, RealVector Ta, RealVector Tb);
extern int Normalize(Real V[], const int dimV, const Real normalizer);
#endif
/* a good practice: end file with a newline */

