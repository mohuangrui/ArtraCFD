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
extern int ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeEigenvaluesAndDecompositionCoefficientAlphaY(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeEigenvaluesAndDecompositionCoefficientAlphaX(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeEigenvaluesAndEigenvectorSpaceLZ(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeEigenvaluesAndEigenvectorSpaceLY(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeEigenvaluesAndEigenvectorSpaceLX(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeEigenvectorSpaceRZ(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeEigenvectorSpaceRY(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeEigenvectorSpaceRX(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
/*
 * Roe average
 *
 * Function
 *      Compute Roe averages.
 */
extern int ComputeRoeAverageZ(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeRoeAverageY(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeRoeAverageX(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
/*
 * Convective fluxes
 *
 * Function
 *      Compute convective fluxes.
 */
extern int ComputeConvectiveFluxZ(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeConvectiveFluxY(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int ComputeConvectiveFluxX(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
extern int sign(const Real x);
extern Real min(const Real x, const Real y);
extern Real max(const Real x, const Real y);
#endif
/* a good practice: end file with a newline */

