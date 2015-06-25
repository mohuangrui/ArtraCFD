/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
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
extern int DecompositionCoefficientAlpha(const int s, Real alpha[],
        const Real deltaU[], const Real Uo[], const Real gamma);
extern int EigenvalueLambda(const int s, Real lambda[], const Real Uo[]);
extern int EigenvectorSpaceLZ(Real L[][DIMU], const Real Uo[], const Real gamma);
extern int EigenvectorSpaceLY(Real L[][DIMU], const Real Uo[], const Real gamma);
extern int EigenvectorSpaceLX(Real L[][DIMU], const Real Uo[], const Real gamma);
extern int EigenvectorSpaceRZ(Real R[][DIMU], const Real Uo[], const Real gamma);
extern int EigenvectorSpaceRY(Real R[][DIMU], const Real Uo[], const Real gamma);
extern int EigenvectorSpaceRX(Real R[][DIMU], const Real Uo[], const Real gamma);
/*
 * Roe average
 *
 * Function
 *      Compute Roe averages.
 */
int ComputeRoeAverage(Real Uo[], const int idx, const int idxh,
        const Real *U, const Real gamma);
/*
 * Convective fluxes
 *
 * Function
 *      Compute convective fluxes.
 */
extern int ConvectiveFluxZ(Real F[], const int idx, const Real *U, const Real gamma);
extern int ConvectiveFluxY(Real F[], const int idx, const Real *U, const Real gamma);
extern int ConvectiveFluxX(Real F[], const int idx, const Real *U, const Real gamma);
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
extern int IndexMath(const int k, const int j, const int i, const Space *);
/*
 * Index geometry
 *
 * Function
 *      Compute the address pointed to current geometry information.
 */
extern Real *IndexGeometry(const int geoID, const Geometry *geometry);
/*
 * Coordinates transformation
 *
 * Function
 *      transform coordinates between node coordinates and general coordinates.
 * Notice
 *      Be cautious with the validity of any calculated index. It's extremely
 *      necessary to adjust the index into the valid interior region or check
 *      validity of the index to avoid index exceed array bound limits and 
 *      mysterious bugs.
 */
extern int ComputeK(const Real z, const Space *);
extern int ComputeJ(const Real y, const Space *);
extern int ComputeI(const Real x, const Space *);
extern int ValidRegionK(const int k, const Partition *);
extern int ValidRegionJ(const int j, const Partition *);
extern int ValidRegionI(const int i, const Partition *);
extern Real ComputeZ(const int k, const Space *);
extern Real ComputeY(const int j, const Space *);
extern Real ComputeX(const int i, const Space *);
/*
 * Common math functions
 */
extern Real MinReal(const Real x, const Real y);
extern Real MaxReal(const Real x, const Real y);
extern int MinInt(const int x, const int y);
extern int MaxInt(const int x, const int y);
extern int Sign(const Real x);
#endif
/* a good practice: end file with a newline */

