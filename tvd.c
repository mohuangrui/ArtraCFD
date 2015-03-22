/****************************************************************************
 * Numeric Scheme for Space Domain                                          *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines the numeric schemes of space domain.                 *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "tvd.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "boundarycondition.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int Lz(Real *U, const Real *Un, const Space *space, 
        const Partition *part, const Flow *flow, const Real dt);
static int Ly(Real *U, const Real *Un, const Space *space, 
        const Partition *part, const Flow *flow, const Real dt);
static int Lx(Real *U, const Real *Un, const Space *space, 
        const Partition *part, const Flow *flow, const Real dt);
static int ComputeReconstructedFluxZ(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeReconstructedFluxY(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeReconstructedFluxX(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int CalculateReconstructedFlux(
        Real Fhat[], const Real F[], const Real Fh[], Real R[][5], const Real Phi[]);
static int ComputeFluxDecompositionCoefficientPhiZ(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeFluxDecompositionCoefficientPhiY(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeFluxDecompositionCoefficientPhiX(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeFunctionGZ(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeFunctionGY(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeFunctionGX(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int CalculateGamma(
        Real gamma[], const Real g[], const Real gh[], const Real alpha[]);
static int CalculateSigma(
        Real sigma[], const Real lambda[], const Real delta[], const Real r);
static int ComputeNumericalDissipationDeltaZ(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Flow *flow);
static int ComputeNumericalDissipationDeltaY(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Flow *flow);
static int ComputeNumericalDissipationDeltaX(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Flow *flow);
static Real Q(const Real z, const Real delta);
static Real minmod(const Real x, const Real y);
static int sign(const Real x);
static Real min(const Real x, const Real y);
static Real max(const Real x, const Real y);
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaY(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaX(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int CalculateAlpha(
        Real alpha[], Real L[][5], const Real deltaU[]);
static int ComputeEigenvaluesAndEigenvectorSpaceLZ(
        Real lambda[], Real L[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeEigenvaluesAndEigenvectorSpaceLY(
        Real lambda[], Real L[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeEigenvaluesAndEigenvectorSpaceLX(
        Real lambda[], Real L[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeEigenvectorSpaceRZ(
        Real R[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeEigenvectorSpaceRY(
        Real R[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeEigenvectorSpaceRX(
        Real R[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeRoeAverageZ(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeRoeAverageY(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeRoeAverageX(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int CalculateRoeAverageUo(
        Real Uo[], const int idx, const int idxh, 
        const Real *U, const Flow *flow);
static int ComputeNonViscousFluxZ(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeNonViscousFluxY(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeNonViscousFluxX(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeViscousFluxGradientZ(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeViscousFluxGradientY(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeViscousFluxGradientX(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeViscousFluxZ(
        Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeViscousFluxY(
        Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeViscousFluxX(
        Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int SpatialDiscretizationAndComputation(Real *U, const Real dt, Real *Uswap, 
        const Space *space, const Particle *particle, 
        const Partition *part, const Flow *flow)
{
    Real *exchanger = U;
    /*
     * When exchange a large bunch of data between two storage space, such as
     * arrays, if there is no new data generation but just data exchange and 
     * update, then the rational way is to exchange the address value that
     * their pointer point to rather than values of data entries.
     */
    Lz(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, particle, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Ly(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, particle, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Lx(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, particle, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Lx(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, particle, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Ly(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, particle, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Lz(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, particle, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */
    /*
     * NOTICE: the exchange operations must be even times, otherwise, the
     * storage space that U originally pointed to will not store the newest
     * updated field data after computation.
     */
    return 0;
}
static int Lz(Real *U, const Real *Un, const Space *space, const Partition *part, const Flow *flow, const Real dt)
{
    Real Fhat[5] = {0.0}; /* reconstructed flux vector */
    Real Fhath[5] = {0.0}; /* reconstructed flux vector at neighbour */
    Real gradG[5] = {0.0}; /* spatial gradient of viscous flux vector */
    int idx = 0; /* linear array index math variable */
    const Real r = dt * space->ddz;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (0 != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                ComputeReconstructedFluxZ(Fhat, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxZ(Fhath, k - 1, j, i, Un, space, flow, dt);
                ComputeViscousFluxGradientZ(gradG, k, j, i, Un, space, flow);
                idx = idx * space->dimU; /* change idx to field variable */
                for (int dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int Ly(Real *U, const Real *Un, const Space *space, const Partition *part, const Flow *flow, const Real dt)
{
    Real Fhat[5] = {0.0}; /* reconstructed flux vector */
    Real Fhath[5] = {0.0}; /* reconstructed flux vector at neighbour */
    Real gradG[5] = {0.0}; /* spatial gradient of viscous flux vector */
    int idx = 0; /* linear array index math variable */
    const Real r = dt * space->ddy;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (0 != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                ComputeReconstructedFluxY(Fhat, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxY(Fhath, k, j - 1, i, Un, space, flow, dt);
                ComputeViscousFluxGradientY(gradG, k, j, i, Un, space, flow);
                idx = idx * space->dimU; /* change idx to field variable */
                for (int dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int Lx(Real *U, const Real *Un, const Space *space, const Partition *part, const Flow *flow, const Real dt)
{
    Real Fhat[5] = {0.0}; /* reconstructed flux vector */
    Real Fhath[5] = {0.0}; /* reconstructed flux vector at neighbour */
    Real gradG[5] = {0.0}; /* spatial gradient of viscous flux vector */
    int idx = 0; /* linear array index math variable */
    const Real r = dt * space->ddx;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (0 != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                ComputeReconstructedFluxX(Fhat, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxX(Fhath, k, j, i - 1, Un, space, flow, dt);
                ComputeViscousFluxGradientX(gradG, k, j, i, Un, space, flow);
                idx = idx * space->dimU; /* change idx to field variable */
                for (int dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int ComputeReconstructedFluxZ(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real F[5] = {0.0}; /* flux at current node */
    Real Fh[5] = {0.0}; /* flux at neighbour */
    Real R[5][5] = {{0.0}}; /* vector space {Rn} */
    Real Phi[5] = {0.0}; /* flux projection or decomposition coefficients on vector space {Rn} */
    ComputeNonViscousFluxZ(F, k, j, i, U, space, flow);
    ComputeNonViscousFluxZ(Fh, k + 1, j, i, U, space, flow);
    ComputeEigenvectorSpaceRZ(R, k, j, i, U, space, flow);
    ComputeFluxDecompositionCoefficientPhiZ(Phi, k, j, i, U, space, flow, dt);
    CalculateReconstructedFlux(Fhat, F, Fh, R, Phi);
    return 0;
}
static int ComputeReconstructedFluxY(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real F[5] = {0.0}; /* flux at current node */
    Real Fh[5] = {0.0}; /* flux at neighbour */
    Real R[5][5] = {{0.0}}; /* vector space {Rn} */
    Real Phi[5] = {0.0}; /* flux projection or decomposition coefficients on vector space {Rn} */
    ComputeNonViscousFluxY(F, k, j, i, U, space, flow);
    ComputeNonViscousFluxY(Fh, k, j + 1, i, U, space, flow);
    ComputeEigenvectorSpaceRY(R, k, j, i, U, space, flow);
    ComputeFluxDecompositionCoefficientPhiY(Phi, k, j, i, U, space, flow, dt);
    CalculateReconstructedFlux(Fhat, F, Fh, R, Phi);
    return 0;
}
static int ComputeReconstructedFluxX(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real F[5] = {0.0}; /* flux at current node */
    Real Fh[5] = {0.0}; /* flux at neighbour */
    Real R[5][5] = {{0.0}}; /* vector space {Rn} */
    Real Phi[5] = {0.0}; /* flux projection or decomposition coefficients on vector space {Rn} */
    ComputeNonViscousFluxX(F, k, j, i, U, space, flow);
    ComputeNonViscousFluxX(Fh, k, j, i + 1, U, space, flow);
    ComputeEigenvectorSpaceRX(R, k, j, i, U, space, flow);
    ComputeFluxDecompositionCoefficientPhiX(Phi, k, j, i, U, space, flow, dt);
    CalculateReconstructedFlux(Fhat, F, Fh, R, Phi);
    return 0;
}
static int CalculateReconstructedFlux(
        Real Fhat[], const Real F[], const Real Fh[], Real R[][5], const Real Phi[])
{
    Real RPhi[5] = {0.0}; /* R x Phi */
    for (int row = 0; row < 5; ++row) {
        RPhi[row] = 0;
        for (int dummy = 0; dummy < 5; ++dummy) {
            RPhi[row] = RPhi[row] + R[row][dummy] * Phi[dummy];
        }
    }
    for (int row = 0; row < 5; ++row) {
        Fhat[row] = 0.5 * (F[row] + Fh[row] + RPhi[row]);
    }
    return 0;
}
static int ComputeFluxDecompositionCoefficientPhiZ(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real g[5] = {0.0}; /* TVD function g at current node */
    Real gh[5] = {0.0}; /* TVD function g at neighbour */
    Real gamma[5] = {0.0}; /* TVD function gamma */
    Real lambda[5] = {0.0}; /* eigenvalues */
    Real alpha[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[5] = {0.0}; /* numerical dissipation */
    ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(lambda, alpha, k, j, i, U, space, flow);
    ComputeFunctionGZ(g, k, j, i, U, space, flow, dt);
    ComputeFunctionGZ(gh, k + 1, j, i, U, space, flow, dt);
    ComputeNumericalDissipationDeltaZ(delta, k, j, i, U, space, flow);
    CalculateGamma(gamma, g, gh, alpha);
    for (int row = 0; row < 5; ++row) {
        Phi[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
    }
    return 0;
}
static int ComputeFluxDecompositionCoefficientPhiY(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real g[5] = {0.0}; /* TVD function g at current node */
    Real gh[5] = {0.0}; /* TVD function g at neighbour */
    Real gamma[5] = {0.0}; /* TVD function gamma */
    Real lambda[5] = {0.0}; /* eigenvalues */
    Real alpha[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[5] = {0.0}; /* numerical dissipation */
    ComputeEigenvaluesAndDecompositionCoefficientAlphaY(lambda, alpha, k, j, i, U, space, flow);
    ComputeFunctionGY(g, k, j, i, U, space, flow, dt);
    ComputeFunctionGY(gh, k, j + 1, i, U, space, flow, dt);
    ComputeNumericalDissipationDeltaY(delta, k, j, i, U, space, flow);
    CalculateGamma(gamma, g, gh, alpha);
    for (int row = 0; row < 5; ++row) {
        Phi[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
    }
    return 0;
}
static int ComputeFluxDecompositionCoefficientPhiX(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real g[5] = {0.0}; /* TVD function g at current node */
    Real gh[5] = {0.0}; /* TVD function g at neighbour */
    Real gamma[5] = {0.0}; /* TVD function gamma */
    Real lambda[5] = {0.0}; /* eigenvalues */
    Real alpha[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[5] = {0.0}; /* numerical dissipation */
    ComputeEigenvaluesAndDecompositionCoefficientAlphaX(lambda, alpha, k, j, i, U, space, flow);
    ComputeFunctionGX(g, k, j, i, U, space, flow, dt);
    ComputeFunctionGX(gh, k, j, i + 1, U, space, flow, dt);
    ComputeNumericalDissipationDeltaX(delta, k, j, i, U, space, flow);
    CalculateGamma(gamma, g, gh, alpha);
    for (int row = 0; row < 5; ++row) {
        Phi[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
    }
    return 0;
}
static int ComputeFunctionGZ(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real lambda[5] = {0.0}; /* eigenvalues */
    Real lambdah[5] = {0.0}; /* eigenvalues at neighbour */
    Real alpha[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real alphah[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[5] = {0.0}; /* numerical dissipation */
    Real deltah[5] = {0.0}; /* numerical dissipation */
    Real sigma[5] = {0.0}; /* TVD function sigma */
    Real sigmah[5] = {0.0}; /* TVD function sigma at neighbour */
    const Real r = dt * space->ddz;
    ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(lambda, alpha, k, j, i, U, space, flow);
    ComputeNumericalDissipationDeltaZ(delta, k, j, i, U, space, flow);
    CalculateSigma(sigma, lambda, delta, r);
    ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(lambdah, alphah, k - 1, j, i, U, space, flow);
    ComputeNumericalDissipationDeltaZ(deltah, k - 1, j, i, U, space, flow);
    CalculateSigma(sigmah, lambdah, deltah, r);
    for (int row = 0; row < 5; ++row) {
        g[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
    }
    return 0;
}
static int ComputeFunctionGY(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real lambda[5] = {0.0}; /* eigenvalues */
    Real lambdah[5] = {0.0}; /* eigenvalues at neighbour */
    Real alpha[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real alphah[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[5] = {0.0}; /* numerical dissipation */
    Real deltah[5] = {0.0}; /* numerical dissipation */
    Real sigma[5] = {0.0}; /* TVD function sigma */
    Real sigmah[5] = {0.0}; /* TVD function sigma at neighbour */
    const Real r = dt * space->ddy;
    ComputeEigenvaluesAndDecompositionCoefficientAlphaY(lambda, alpha, k, j, i, U, space, flow);
    ComputeNumericalDissipationDeltaY(delta, k, j, i, U, space, flow);
    CalculateSigma(sigma, lambda, delta, r);
    ComputeEigenvaluesAndDecompositionCoefficientAlphaY(lambdah, alphah, k, j - 1, i, U, space, flow);
    ComputeNumericalDissipationDeltaY(deltah, k, j - 1, i, U, space, flow);
    CalculateSigma(sigmah, lambdah, deltah, r);
    for (int row = 0; row < 5; ++row) {
        g[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
    }
    return 0;
}
static int ComputeFunctionGX(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real lambda[5] = {0.0}; /* eigenvalues */
    Real lambdah[5] = {0.0}; /* eigenvalues at neighbour */
    Real alpha[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real alphah[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[5] = {0.0}; /* numerical dissipation */
    Real deltah[5] = {0.0}; /* numerical dissipation */
    Real sigma[5] = {0.0}; /* TVD function sigma */
    Real sigmah[5] = {0.0}; /* TVD function sigma at neighbour */
    const Real r = dt * space->ddx;
    ComputeEigenvaluesAndDecompositionCoefficientAlphaX(lambda, alpha, k, j, i, U, space, flow);
    ComputeNumericalDissipationDeltaX(delta, k, j, i, U, space, flow);
    CalculateSigma(sigma, lambda, delta, r);
    ComputeEigenvaluesAndDecompositionCoefficientAlphaX(lambdah, alphah, k, j, i - 1, U, space, flow);
    ComputeNumericalDissipationDeltaX(deltah, k, j, i - 1, U, space, flow);
    CalculateSigma(sigmah, lambdah, deltah, r);
    for (int row = 0; row < 5; ++row) {
        g[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
    }
    return 0;
}
static int CalculateGamma(
        Real gamma[], const Real g[], const Real gh[], const Real alpha[])
{
    for (int row = 0; row < 5; ++row) {
        if (0 != alpha[row]) {
            gamma[row] = (gh[row] - g[row]) / alpha[row];
        }
        else {
            gamma[row] = 0;
        }
    }
    return 0;
}
static int CalculateSigma(
        Real sigma[], const Real lambda[], const Real delta[], const Real r)
{
    for (int row = 0; row < 5; ++row) {
        sigma[row] = 0.5 * (Q(lambda[row], delta[row]) - r * lambda[row] * lambda[row]);
    }
    return 0;
}
static int ComputeNumericalDissipationDeltaZ(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    /* numerical dissipation in [0.05, 0.25], 0.125 is recommended */
    ComputeRoeAverageZ(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    for (int row = 0; row < 5; ++row) {
        delta[row] = flow->delta * (fabs(u) + fabs(v) + fabs(w) + c); 
    }
    return 0;
}
static int ComputeNumericalDissipationDeltaY(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    /* numerical dissipation in [0.05, 0.25], 0.125 is recommended */
    ComputeRoeAverageY(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    for (int row = 0; row < 5; ++row) {
        delta[row] = flow->delta * (fabs(u) + fabs(v) + fabs(w) + c); 
    }
    return 0;
}
static int ComputeNumericalDissipationDeltaX(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    /* numerical dissipation in [0.05, 0.25], 0.125 is recommended */
    ComputeRoeAverageX(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    for (int row = 0; row < 5; ++row) {
        delta[row] = flow->delta * (fabs(u) + fabs(v) + fabs(w) + c); 
    }
    return 0;
}
static Real Q(const Real z, const Real delta)
{
    if (fabs(z) >= delta) {
        return fabs(z);
    }
    return (0.5 * (z * z / delta + delta));
}
static Real minmod(const Real x, const Real y)
{
    return (sign(x) * max(0, min(fabs(x), y * sign(x))));
}
static int sign(const Real x)
{
    if (0 < x) {
        return 1;
    }
    if (0 > x) {
        return -1;
    }
    return 0;
}
static Real min(const Real x, const Real y)
{
    if (x < y) {
        return x;
    }
    return y;
}
static Real max(const Real x, const Real y)
{
    if (x > y) {
        return x;
    }
    return y;
}
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxh = IndexMath(k + 1, j, i, space) * space->dimU;;
    Real L[5][5] = {{0.0}}; /* store left eigenvectors */
    const Real deltaU[5] = {
        U[idxh+0] - U[idx+0],
        U[idxh+1] - U[idx+1],
        U[idxh+2] - U[idx+2],
        U[idxh+3] - U[idx+3],
        U[idxh+4] - U[idx+4]};
    ComputeEigenvaluesAndEigenvectorSpaceLZ(lambda, L, k, j, i, U, space, flow);
    CalculateAlpha(alpha, L, deltaU);
    return 0;
}
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaY(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxh = IndexMath(k, j + 1, i, space) * space->dimU;;
    Real L[5][5] = {{0.0}}; /* store left eigenvectors */
    const Real deltaU[5] = {
        U[idxh+0] - U[idx+0],
        U[idxh+1] - U[idx+1],
        U[idxh+2] - U[idx+2],
        U[idxh+3] - U[idx+3],
        U[idxh+4] - U[idx+4]};
    ComputeEigenvaluesAndEigenvectorSpaceLY(lambda, L, k, j, i, U, space, flow);
    CalculateAlpha(alpha, L, deltaU);
    return 0;
}
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaX(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxh = IndexMath(k, j, i + 1, space) * space->dimU;;
    Real L[5][5] = {{0.0}}; /* store left eigenvectors */
    const Real deltaU[5] = {
        U[idxh+0] - U[idx+0],
        U[idxh+1] - U[idx+1],
        U[idxh+2] - U[idx+2],
        U[idxh+3] - U[idx+3],
        U[idxh+4] - U[idx+4]};
    ComputeEigenvaluesAndEigenvectorSpaceLX(lambda, L, k, j, i, U, space, flow);
    CalculateAlpha(alpha, L, deltaU);
    return 0;
}
static int CalculateAlpha(
        Real alpha[], Real L[][5], const Real deltaU[])
{
    for (int row = 0; row < 5; ++row) {
        alpha[row] = 0;
        for (int dummy = 0; dummy < 5; ++dummy) {
            alpha[row] = alpha[row] + L[row][dummy] * deltaU[dummy];
        }
    }
    return 0;
}
static int ComputeEigenvaluesAndEigenvectorSpaceLZ(
        Real lambda[], Real L[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageZ(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (flow->gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    lambda[0] = w - c; lambda[1] = w; lambda[2] = w; lambda[3] = w; lambda[4] = w + c;
    L[0][0] = b * q + d * w;   L[0][1] = -b * u;             L[0][2] = -b * v;             L[0][3] = -b * w - d;     L[0][4] = b;
    L[1][0] = -2 * b * q * u;  L[1][1] = 2 * b * u * u + 1;  L[1][2] = 2 * b * v * u;      L[1][3] = 2 * b * w * u;  L[1][4] = -2 * b * u;
    L[2][0] = -2 * b * q * v;  L[2][1] = 2 * b * v * u;      L[2][2] = 2 * b * v * v + 1;  L[2][3] = 2 * b * w * v;  L[2][4] = -2 * b * v;
    L[3][0] = -2 * b * q + 1;  L[3][1] = 2 * b * u;          L[3][2] = 2 * b * v;          L[3][3] = 2 * b * w;      L[3][4] = -2 * b;
    L[4][0] = b * q - d * w;   L[4][1] = -b * u;             L[4][2] = -b * v;             L[4][3] = -b * w + d;     L[4][4] = b;
    return 0;
}
static int ComputeEigenvaluesAndEigenvectorSpaceLY(
        Real lambda[], Real L[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageY(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (flow->gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    lambda[0] = v - c; lambda[1] = v; lambda[2] = v; lambda[3] = v; lambda[4] = v + c;
    L[0][0] = b * q + d * v;    L[0][1] = -b * u;             L[0][2] = -b * v - d;     L[0][3] = -b * w;             L[0][4] = b;
    L[1][0] = -2 * b * q * u;   L[1][1] = 2 * b * u * u + 1;  L[1][2] = 2 * b * v * u;  L[1][3] = 2 * b * w * u;      L[1][4] = -2 * b * u;
    L[2][0] = -2 * b * q + 1;   L[2][1] = 2 * b * u;          L[2][2] = 2 * b * v;      L[2][3] = 2 * b * w;          L[2][4] = -2 * b;
    L[3][0] = -2 * b * q * w;   L[3][1] = 2 * b * w * u;      L[3][2] = 2 * b * w * v;  L[3][3] = 2 * b * w * w + 1;  L[3][4] = -2 * b * w;
    L[4][0] = b * q - d * v;    L[4][1] = -b * u;             L[4][2] = -b * v + d;     L[4][3] = -b * w;             L[4][4] = b;
    return 0;
}
static int ComputeEigenvaluesAndEigenvectorSpaceLX(
        Real lambda[], Real L[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageX(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (flow->gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    lambda[0] = u - c; lambda[1] = u; lambda[2] = u; lambda[3] = u; lambda[4] = u + c;
    L[0][0] = b * q + d * u;   L[0][1] = -b * u - d;     L[0][2] = -b * v;             L[0][3] = -b * w;             L[0][4] = b;
    L[1][0] = -2 * b * q + 1;  L[1][1] = 2 * b * u;      L[1][2] = 2 * b * v;          L[1][3] = 2 * b * w;          L[1][4] = -2 * b;
    L[2][0] = -2 * b * q * v;  L[2][1] = 2 * b * v * u;  L[2][2] = 2 * b * v * v + 1;  L[2][3] = 2 * b * w * v;      L[2][4] = -2 * b * v;
    L[3][0] = -2 * b * q * w;  L[3][1] = 2 * b * w * u;  L[3][2] = 2 * b * w * v;      L[3][3] = 2 * b * w * w + 1;  L[3][4] = -2 * b * w;
    L[4][0] = b * q - d * u;   L[4][1] = -b * u + d;     L[4][2] = -b * v;             L[4][3] = -b * w;             L[4][4] = b;
    return 0;
}
static int ComputeEigenvectorSpaceRZ(
        Real R[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageZ(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real hT = Uo[4];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    R[0][0] = 1;           R[0][1] = 0;  R[0][2] = 0;  R[0][3] = 1;          R[0][4] = 1;
    R[1][0] = u;           R[1][1] = 1;  R[1][2] = 0;  R[1][3] = 0;          R[1][4] = u;
    R[2][0] = v;           R[2][1] = 0;  R[2][2] = 1;  R[2][3] = 0;          R[2][4] = v;
    R[3][0] = w - c;       R[3][1] = 0;  R[3][2] = 0;  R[3][3] = w;          R[3][4] = w + c;
    R[4][0] = hT - w * c;  R[4][1] = u;  R[4][2] = v;  R[4][3] = w * w - q;  R[4][4] = hT + w * c;
    return 0;
}
static int ComputeEigenvectorSpaceRY(
        Real R[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageY(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real hT = Uo[4];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    R[0][0] = 1;           R[0][1] = 0;  R[0][2] = 1;          R[0][3] = 0;  R[0][4] = 1;
    R[1][0] = u;           R[1][1] = 1;  R[1][2] = 0;          R[1][3] = 0;  R[1][4] = u;
    R[2][0] = v - c;       R[2][1] = 0;  R[2][2] = v;          R[2][3] = 0;  R[2][4] = v + c;
    R[3][0] = w;           R[3][1] = 0;  R[3][2] = 0;          R[3][3] = 1;  R[3][4] = w;
    R[4][0] = hT - v * c;  R[4][1] = u;  R[4][2] = v * v - q;  R[4][3] = w;  R[4][4] = hT + v * c;
    return 0;
}
static int ComputeEigenvectorSpaceRX(
        Real R[][5], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageX(Uo, k, j, i, U, space, flow);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real hT = Uo[4];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    R[0][0] = 1;           R[0][1] = 1;          R[0][2] = 0;  R[0][3] = 0;  R[0][4] = 1;
    R[1][0] = u - c;       R[1][1] = u;          R[1][2] = 0;  R[1][3] = 0;  R[1][4] = u + c;
    R[2][0] = v;           R[2][1] = 0;          R[2][2] = 1;  R[2][3] = 0;  R[2][4] = v;
    R[3][0] = w;           R[3][1] = 0;          R[3][2] = 0;  R[3][3] = 1;  R[3][4] = w;
    R[4][0] = hT - u * c;  R[4][1] = u * u - q;  R[4][2] = v;  R[4][3] = w;  R[4][4] = hT + u * c;
    return 0;
}
static int ComputeRoeAverageZ(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxh = IndexMath(k + 1, j, i, space) * space->dimU;;
    CalculateRoeAverageUo(Uo, idx, idxh, U, flow);
    return 0;
}
static int ComputeRoeAverageY(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxh = IndexMath(k, j + 1, i, space) * space->dimU;;
    CalculateRoeAverageUo(Uo, idx, idxh, U, flow);
    return 0;
}
static int ComputeRoeAverageX(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxh = IndexMath(k, j, i + 1, space) * space->dimU;;
    CalculateRoeAverageUo(Uo, idx, idxh, U, flow);
    return 0;
}
static int CalculateRoeAverageUo(
        Real Uo[], const int idx, const int idxh, 
        const Real *U, const Flow *flow)
{
    const Real rho = U[idx+0];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real hT = flow->gamma * U[idx+4] / rho - 0.5 * (u * u + v * v + w * w) * (flow->gamma - 1.0);
    const Real rho_h = U[idxh+0];
    const Real u_h = U[idxh+1] / rho_h;
    const Real v_h = U[idxh+2] / rho_h;
    const Real w_h = U[idxh+3] / rho_h;
    const Real hT_h = flow->gamma * U[idxh+4] / rho_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h) * (flow->gamma - 1.0);
    const Real D = sqrt(rho_h / rho);
    //Uo[0] = rho * (0.5 * (1.0 + D)) * (0.5 * (1.0 + D)); /* rho average, not required */
    Uo[1] = (u + D * u_h) / (1.0 + D); /* u average */
    Uo[2] = (v + D * v_h) / (1.0 + D); /* v average */
    Uo[3] = (w + D * w_h) / (1.0 + D); /* w average */
    Uo[4] = (hT + D * hT_h) / (1.0 + D); /* hT average */
    Uo[5] = sqrt((flow->gamma - 1.0) * (Uo[4] - 0.5 * (Uo[1] * Uo[1] + Uo[2] * Uo[2] + Uo[3] * Uo[3]))); /* the speed of sound */
    return 0;
}
static int ComputeNonViscousFluxZ(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const Real rho = U[idx+0];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = (flow->gamma - 1.0) * rho * (eT - 0.5 * (u * u + v * v + w * w));
    F[0] = rho * w;
    F[1] = rho * w * u;
    F[2] = rho * w * v;
    F[3] = rho * w * w + p;
    F[4] = (rho * eT + p) * w;
    return 0;
}
static int ComputeNonViscousFluxY(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const Real rho = U[idx+0];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = (flow->gamma - 1.0) * rho * (eT - 0.5 * (u * u + v * v + w * w));
    F[0] = rho * v;
    F[1] = rho * v * u;
    F[2] = rho * v * v + p;
    F[3] = rho * v * w;
    F[4] = (rho * eT + p) * v;
    return 0;
}
static int ComputeNonViscousFluxX(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const Real rho = U[idx+0];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = (flow->gamma - 1.0) * rho * (eT - 0.5 * (u * u + v * v + w * w));
    F[0] = rho * u;
    F[1] = rho * u * u + p;
    F[2] = rho * u * v;
    F[3] = rho * u * w;
    F[4] = (rho * eT + p) * u;
    return 0;
}
/*
 * Generally the viscous terms will only be discretized by central
 * difference scheme. However, for outermost layer of interior nodes, the
 * central scheme of them require the viscous flux vector at boundaries, 
 * because the corners of current computational domain haven't set with any
 * values, the viscous flux vector at boundaries can not be interpolated
 * with central scheme, therefore, need to change to forward or backward. 
 * This situation is the same to interior ghost nodes the central
 * difference scheme can't be applied because of lacking stencil. Thus,
 * they also need to be identified and modified.
 */
static int ComputeViscousFluxGradientZ(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real hG[5] = {0.0}; /* viscous flux vector */
    Real Gh[5] = {0.0}; /* viscous flux vector */
    Real h = 0; /* reciprocal of differencing distance */
    const int offset = space->nodeFlagOffset;
    /* default is central scheme */
    int hl = k - 1;
    int hr = k + 1;
    if (space->ng == hl) { /* if left at boundary */
        ++hl;
    }
    if ((space->nz + space->ng - 1) == hr) { /* if right at boundary */
        --hr;
    }
    /* check ghost */
    const int idxl = IndexMath(hl, j, i, space);
    const int idxr = IndexMath(hr, j, i, space);
    if (offset <= space->nodeFlag[idxl]) {
        ++hl;
    }
    if (offset <= space->nodeFlag[idxr]) {
        --hr;
    }
    if (0 != (hr - hl)) { /* only do calculation when needed */
        h = space->ddz / (hr - hl);
        ComputeViscousFluxZ(hG, hl, j, i, U, space, flow);
        ComputeViscousFluxZ(Gh, hr, j, i, U, space, flow);
    }
    for (int row = 0; row < 5; ++row) {
        gradG[row] = h * (Gh[row] - hG[row]);
    }
    return 0;
}
static int ComputeViscousFluxGradientY(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real hG[5] = {0.0}; /* viscous flux vector */
    Real Gh[5] = {0.0}; /* viscous flux vector */
    Real h = 0; /* reciprocal of differencing distance */
    const int offset = space->nodeFlagOffset;
    /* default is central scheme */
    int hl = j - 1;
    int hr = j + 1;
    if (space->ng == hl) { /* if left at boundary */
        ++hl;
    }
    if ((space->ny + space->ng - 1) == hr) { /* if right at boundary */
        --hr;
    }
    /* check ghost */
    const int idxl = IndexMath(k, hl, i, space);
    const int idxr = IndexMath(k, hr, i, space);
    if (offset <= space->nodeFlag[idxl]) {
        ++hl;
    }
    if (offset <= space->nodeFlag[idxr]) {
        --hr;
    }
    if (0 != (hr - hl)) { /* only do calculation when needed */
        h = space->ddy / (hr - hl);
        ComputeViscousFluxY(hG, k, hl, i, U, space, flow);
        ComputeViscousFluxY(Gh, k, hr, i, U, space, flow);
    }
    for (int row = 0; row < 5; ++row) {
        gradG[row] = h * (Gh[row] - hG[row]);
    }
    return 0;
}
static int ComputeViscousFluxGradientX(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real hG[5] = {0.0}; /* viscous flux vector */
    Real Gh[5] = {0.0}; /* viscous flux vector */
    Real h = 0; /* reciprocal of differencing distance */
    const int offset = space->nodeFlagOffset;
    /* default is central scheme */
    int hl = i - 1;
    int hr = i + 1;
    if (space->ng == hl) { /* if left at boundary */
        ++hl;
    }
    if ((space->nx + space->ng - 1) == hr) { /* if right at boundary */
        --hr;
    }
    /* check ghost */
    const int idxl = IndexMath(k, j, hl, space);
    const int idxr = IndexMath(k, j, hr, space);
    if (offset <= space->nodeFlag[idxl]) {
        ++hl;
    }
    if (offset <= space->nodeFlag[idxr]) {
        --hr;
    }
    if (0 != (hr - hl)) { /* only do calculation when needed */
        h = space->ddx / (hr - hl);
        ComputeViscousFluxX(hG, k, j, hl, U, space, flow);
        ComputeViscousFluxX(Gh, k, j, hr, U, space, flow);
    }
    for (int row = 0; row < 5; ++row) {
        gradG[row] = h * (Gh[row] - hG[row]);
    }
    return 0;
}
static int ComputeViscousFluxZ(
        Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxW = IndexMath(k, j, i - 1, space) * space->dimU;
    const int idxE = IndexMath(k, j, i + 1, space) * space->dimU;
    const int idxS = IndexMath(k, j - 1, i, space) * space->dimU;
    const int idxN = IndexMath(k, j + 1, i, space) * space->dimU;
    const int idxF = IndexMath(k - 1, j, i, space) * space->dimU;
    const int idxB = IndexMath(k + 1, j, i, space) * space->dimU;

    /* calculate derivatives in z direction */
    const Real rhoB = U[idxB+0];
    const Real uB = U[idxB+1] / rhoB;
    const Real vB = U[idxB+2] / rhoB;
    const Real wB = U[idxB+3] / rhoB;
    const Real eTB = U[idxB+4] / rhoB;
    const Real TB = (eTB - 0.5 * (uB * uB + vB * vB + wB * wB)) / flow->cv;

    const Real rhoF = U[idxF+0];
    const Real uF = U[idxF+1] / rhoF;
    const Real vF = U[idxF+2] / rhoF;
    const Real wF = U[idxF+3] / rhoF;
    const Real eTF = U[idxF+4] / rhoF;
    const Real TF = (eTF - 0.5 * (uF * uF + vF * vF + wF * wF)) / flow->cv;

    const Real du_dz = (uB - uF) * (0.5 * space->ddz);
    const Real dv_dz = (vB - vF) * (0.5 * space->ddz);
    const Real dw_dz = (wB - wF) * (0.5 * space->ddz);
    const Real dT_dz = (TB - TF) * (0.5 * space->ddz);

    /* calculate derivatives in y direction */
    const Real vN = U[idxN+2] / U[idxN+0];
    const Real wN = U[idxN+3] / U[idxN+0];
    const Real vS = U[idxS+2] / U[idxS+0];
    const Real wS = U[idxS+3] / U[idxS+0];
    const Real dv_dy = (vN - vS) * (0.5 * space->ddy);
    const Real dw_dy = (wN - wS) * (0.5 * space->ddy);

    /* calculate derivatives in x direction */
    const Real uE = U[idxE+1] / U[idxE+0];
    const Real wE = U[idxE+3] / U[idxE+0];
    const Real uW = U[idxW+1] / U[idxW+0];
    const Real wW = U[idxW+3] / U[idxW+0];
    const Real du_dx = (uE - uW) * (0.5 * space->ddx);
    const Real dw_dx = (wE - wW) * (0.5 * space->ddx);

    /* the primitive variables in current point */
    const Real rho = U[idx+0];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

    /* Calculate dynamic viscosity and heat conductivity */
    const Real mu = flow->refMu * 1.45e-6 * (pow(T * flow->refTemperature, 1.5) / (T * flow->refTemperature + 110));
    const Real heatK = flow->gamma * flow->cv * mu / flow->refPr;
    const Real divV = du_dx + dv_dy + dw_dz;

    G[0] = 0;
    G[1] = mu * (dw_dx + du_dz);
    G[2] = mu * (dw_dy + dv_dz);
    G[3] = mu * (2.0 * dw_dz - (2.0/3.0) * divV);
    G[4] = heatK * dT_dz + u * G[1] + v * G[2] + w * G[3];
    return 0;
}
static int ComputeViscousFluxY(
        Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxW = IndexMath(k, j, i - 1, space) * space->dimU;
    const int idxE = IndexMath(k, j, i + 1, space) * space->dimU;
    const int idxS = IndexMath(k, j - 1, i, space) * space->dimU;
    const int idxN = IndexMath(k, j + 1, i, space) * space->dimU;
    const int idxF = IndexMath(k - 1, j, i, space) * space->dimU;
    const int idxB = IndexMath(k + 1, j, i, space) * space->dimU;

    /* calculate derivatives in z direction */
    const Real vB = U[idxB+2] / U[idxB+0];
    const Real wB = U[idxB+3] / U[idxB+0];
    const Real vF = U[idxF+2] / U[idxF+0];
    const Real wF = U[idxF+3] / U[idxF+0];
    const Real dv_dz = (vB - vF) * (0.5 * space->ddz);
    const Real dw_dz = (wB - wF) * (0.5 * space->ddz);

    /* calculate derivatives in y direction */
    const Real rhoN = U[idxN+0];
    const Real uN = U[idxN+1] / rhoN;
    const Real vN = U[idxN+2] / rhoN;
    const Real wN = U[idxN+3] / rhoN;
    const Real eTN = U[idxN+4] / rhoN;
    const Real TN = (eTN - 0.5 * (uN * uN + vN * vN + wN * wN)) / flow->cv;

    const Real rhoS = U[idxS+0];
    const Real uS = U[idxS+1] / rhoS;
    const Real vS = U[idxS+2] / rhoS;
    const Real wS = U[idxS+3] / rhoS;
    const Real eTS = U[idxS+4] / rhoS;
    const Real TS = (eTS - 0.5 * (uS * uS + vS * vS + wS * wS)) / flow->cv;

    const Real du_dy = (uN - uS) * (0.5 * space->ddy);
    const Real dv_dy = (vN - vS) * (0.5 * space->ddy);
    const Real dw_dy = (wN - wS) * (0.5 * space->ddy);
    const Real dT_dy = (TN - TS) * (0.5 * space->ddy);

    /* calculate derivatives in x direction */
    const Real uE = U[idxE+1] / U[idxE+0];
    const Real vE = U[idxE+2] / U[idxE+0];
    const Real uW = U[idxW+1] / U[idxW+0];
    const Real vW = U[idxW+2] / U[idxW+0];
    const Real du_dx = (uE - uW) * (0.5 * space->ddx);
    const Real dv_dx = (vE - vW) * (0.5 * space->ddx);

    /* the primitive variables in current point */
    const Real rho = U[idx+0];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

    /* Calculate dynamic viscosity and heat conductivity */
    const Real mu = flow->refMu * 1.45e-6 * (pow(T * flow->refTemperature, 1.5) / (T * flow->refTemperature + 110));
    const Real heatK = flow->gamma * flow->cv * mu / flow->refPr;
    const Real divV = du_dx + dv_dy + dw_dz;

    G[0] = 0;
    G[1] = mu * (dv_dx + du_dy);
    G[2] = mu * (2.0 * dv_dy - (2.0/3.0) * divV);
    G[3] = mu * (dv_dz + dw_dy);
    G[4] = heatK * dT_dy + u * G[1] + v * G[2] + w * G[3];
    return 0;
}
static int ComputeViscousFluxX(
        Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = IndexMath(k, j, i, space) * space->dimU;
    const int idxW = IndexMath(k, j, i - 1, space) * space->dimU;
    const int idxE = IndexMath(k, j, i + 1, space) * space->dimU;
    const int idxS = IndexMath(k, j - 1, i, space) * space->dimU;
    const int idxN = IndexMath(k, j + 1, i, space) * space->dimU;
    const int idxF = IndexMath(k - 1, j, i, space) * space->dimU;
    const int idxB = IndexMath(k + 1, j, i, space) * space->dimU;

    /* calculate derivatives in z direction */
    const Real uB = U[idxB+1] / U[idxB+0];
    const Real uF = U[idxF+1] / U[idxF+0];
    const Real wB = U[idxB+3] / U[idxB+0];
    const Real wF = U[idxF+3] / U[idxF+0];
    const Real du_dz = (uB - uF) * (0.5 * space->ddz);
    const Real dw_dz = (wB - wF) * (0.5 * space->ddz);

    /* calculate derivatives in y direction */
    const Real uN = U[idxN+1] / U[idxN+0];
    const Real uS = U[idxS+1] / U[idxS+0];
    const Real vN = U[idxN+2] / U[idxN+0];
    const Real vS = U[idxS+2] / U[idxS+0];
    const Real du_dy = (uN - uS) * (0.5 * space->ddy);
    const Real dv_dy = (vN - vS) * (0.5 * space->ddy);

    /* calculate derivatives in x direction */
    const Real rhoE = U[idxE+0];
    const Real uE = U[idxE+1] / rhoE;
    const Real vE = U[idxE+2] / rhoE;
    const Real wE = U[idxE+3] / rhoE;
    const Real eTE = U[idxE+4] / rhoE;
    const Real TE = (eTE - 0.5 * (uE * uE + vE * vE + wE * wE)) / flow->cv;

    const Real rhoW = U[idxW+0];
    const Real uW = U[idxW+1] / rhoW;
    const Real vW = U[idxW+2] / rhoW;
    const Real wW = U[idxW+3] / rhoW;
    const Real eTW = U[idxW+4] / rhoW;
    const Real TW = (eTW - 0.5 * (uW * uW + vW * vW + wW * wW)) / flow->cv;

    const Real du_dx = (uE - uW) * (0.5 * space->ddx);
    const Real dv_dx = (vE - vW) * (0.5 * space->ddx);
    const Real dw_dx = (wE - wW) * (0.5 * space->ddx);
    const Real dT_dx = (TE - TW) * (0.5 * space->ddx);

    /* the primitive variables in current point */
    const Real rho = U[idx+0];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

    /* Calculate dynamic viscosity and heat conductivity */
    const Real mu = flow->refMu * 1.45e-6 * (pow(T * flow->refTemperature, 1.5) / (T * flow->refTemperature + 110));
    const Real heatK = flow->gamma * flow->cv * mu / flow->refPr;
    const Real divV = du_dx + dv_dy + dw_dz;

    G[0] = 0;
    G[1] = mu * (2.0 * du_dx - (2.0/3.0) * divV);
    G[2] = mu * (du_dy + dv_dx);
    G[3] = mu * (du_dz + dw_dx);
    G[4] = heatK * dT_dx + u * G[1] + v * G[2] + w * G[3];
    return 0;
}
/* a good practice: end file with a newline */

