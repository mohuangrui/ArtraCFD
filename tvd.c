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
static int ComputeReconstructedFluxTVD(
        Real Fhatz[], Real Fhaty[], Real Fhatx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int CalculateReconstructedFlux(
        Real Fhat[], const Real F[], const Real Fh[], Real R[][5], const Real Phi[]);
static int ComputeEigenvectorSpaceR(
        Real Rz[][5], Real Ry[][5], Real Rx[][5], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeFluxDecompositionCoefficientPhi(
        Real Phiz[], Real Phiy[], Real Phix[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeEigenvaluesAndDecompositionCoefficientAlpha(
        Real lambdaz[], Real lambday[], Real lambdax[], 
        Real alphaz[], Real alphay[], Real alphax[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int CalculateAlpha(
        Real alpha[], Real L[][5], const Real deltaU[]);
static int ComputeEigenvaluesAndEigenvectorSpaceL(
        Real lambdaz[], Real lambday[], Real lambdax[],
        Real Lz[][5], Real Ly[][5], Real Lx[][5], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeFunctionG(
        Real gz[], Real gy[], Real gx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int CalculateGamma(
        Real gamma[], const Real g[], const Real gh[], const Real alpha[]);
static int CalculateSigma(
        Real sigma[], const Real lambda[], const Real delta[], const Real r);
static int ComputeRoeAverage(Real Uoz[], Real Uoy[], Real Uox[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int CalculateRoeAverageUo(
        Real Uo[], const int idx, const int idxh, 
        const Real *U, const Flow *flow);
static int ComputeNonViscousFlux(Real Fz[], Real Fy[], Real Fx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeViscousFluxGradient(Real gradGz[], Real gradGy[], Real gradGx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeViscousFlux(Real Gz[], Real Gy[], Real Gx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow);
static int ComputeNumericalDissipationDelta(
        Real deltaz[], Real deltay[], Real deltax[], 
        const int k, const int j, const int i,
        const Real *U, const Space *space, const Flow *flow);
static Real Q(const Real z, const Real delta);
static Real minmod(const Real x, const Real y);
static int sign(const Real x);
static Real min(const Real x, const Real y);
static Real max(const Real x, const Real y);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int SpatialDiscretizationAndComputation(Real *U, Real *Un, 
        const Space *space, const Particle *particle, 
        const Partition *part, const Flow *flow, const Real dt)
{
    Real *exchanger = Un;
    /*
     * When exchange a large bunch of data between two storage space, such as
     * arrays, if there is no new data generation but just data exchange and 
     * update, then the rational way is to exchange the address value that
     * their pointer point to rather than values of data entries.
     */
    Lz(U, Un, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(U, space, particle, part, flow);
    exchanger = Un; /* preserve the address of Un */
    Un = U; /* update flow field */
    U = exchanger; /* regain the used space as new space */

    Ly(U, Un, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(U, space, particle, part, flow);
    exchanger = Un; /* preserve the address of Un */
    Un = U; /* update flow field */
    U = exchanger; /* regain the used space as new space */

    Lx(U, Un, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(U, space, particle, part, flow);
    exchanger = Un; /* preserve the address of Un */
    Un = U; /* update flow field */
    U = exchanger; /* regain the used space as new space */

    Lx(U, Un, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(U, space, particle, part, flow);
    exchanger = Un; /* preserve the address of Un */
    Un = U; /* update flow field */
    U = exchanger; /* regain the used space as new space */

    Ly(U, Un, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(U, space, particle, part, flow);
    exchanger = Un; /* preserve the address of Un */
    Un = U; /* update flow field */
    U = exchanger; /* regain the used space as new space */

    Lz(U, Un, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(U, space, particle, part, flow);
    exchanger = Un; /* preserve the address of Un */
    Un = U; /* update flow field */
    U = exchanger; /* regain the used space as new space */
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
                idx = ((k * space->jMax + j) * space->iMax + i);
                if (0 != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                idx = idx * 5; /* change idx to field variable */
                ComputeReconstructedFluxTVD(Fhat, NULL, NULL, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxTVD(Fhath, NULL, NULL, k - 1, j, i, Un, space, flow, dt);
                ComputeViscousFluxGradient(gradG, NULL, NULL, k, j, i, Un, space, flow);
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
                idx = ((k * space->jMax + j) * space->iMax + i);
                if (0 != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                idx = idx * 5; /* change idx to field variable */
                ComputeReconstructedFluxTVD(NULL, Fhat, NULL, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxTVD(NULL, Fhath, NULL, k, j - 1, i, Un, space, flow, dt);
                ComputeViscousFluxGradient(NULL, gradG, NULL, k, j, i, Un, space, flow);
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
                idx = ((k * space->jMax + j) * space->iMax + i);
                if (0 != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                idx = idx * 5; /* change idx to field variable */
                ComputeReconstructedFluxTVD(NULL, NULL, Fhat, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxTVD(NULL, NULL, Fhath, k, j, i - 1, Un, space, flow, dt);
                ComputeViscousFluxGradient(NULL, NULL, gradG, k, j, i, Un, space, flow);
                for (int dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int ComputeReconstructedFluxTVD(
        Real Fhatz[], Real Fhaty[], Real Fhatx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real F[5] = {0.0}; /* flux at current node */
    Real Fh[5] = {0.0}; /* flux at neighbour */
    Real R[5][5] = {{0.0}}; /* vector space {Rn} */
    Real Phi[5] = {0.0}; /* flux projection or decomposition coefficients on vector space {Rn} */
    if (NULL != Fhatz) {
        ComputeNonViscousFlux(F, NULL, NULL, k, j, i, U, space, flow);
        ComputeNonViscousFlux(Fh, NULL, NULL, k + 1, j, i, U, space, flow);
        ComputeEigenvectorSpaceR(R, NULL, NULL, k, j, i, U, space, flow);
        ComputeFluxDecompositionCoefficientPhi(Phi, NULL, NULL, k, j, i, U, space, flow, dt);
        CalculateReconstructedFlux(Fhatz, F, Fh, R, Phi);
    }
    if (NULL != Fhaty) {
        ComputeNonViscousFlux(NULL, F, NULL, k, j, i, U, space, flow);
        ComputeNonViscousFlux(NULL, Fh, NULL, k, j + 1, i, U, space, flow);
        ComputeEigenvectorSpaceR(NULL, R, NULL, k, j, i, U, space, flow);
        ComputeFluxDecompositionCoefficientPhi(NULL, Phi, NULL, k, j, i, U, space, flow, dt);
        CalculateReconstructedFlux(Fhaty, F, Fh, R, Phi);
    }
    if (NULL != Fhatx) {
        ComputeNonViscousFlux(NULL, NULL, F, k, j, i, U, space, flow);
        ComputeNonViscousFlux(NULL, NULL, Fh, k, j, i + 1, U, space, flow);
        ComputeEigenvectorSpaceR(NULL, NULL, R, k, j, i, U, space, flow);
        ComputeFluxDecompositionCoefficientPhi(NULL, NULL, Phi, k, j, i, U, space, flow, dt);
        CalculateReconstructedFlux(Fhatx, F, Fh, R, Phi);
    }
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
static int ComputeFluxDecompositionCoefficientPhi(
        Real Phiz[], Real Phiy[], Real Phix[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt)
{
    Real g[5] = {0.0}; /* TVD function g at current node */
    Real gh[5] = {0.0}; /* TVD function g at neighbour */
    Real gamma[5] = {0.0}; /* TVD function gamma */
    Real lambda[5] = {0.0}; /* eigenvalues */
    Real alpha[5] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[5] = {0.0}; /* numerical dissipation */
    if (NULL != Phiz) {
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(lambda, NULL, NULL, alpha, NULL, NULL, k, j, i, U, space, flow);
        ComputeFunctionG(g, NULL, NULL, k, j, i, U, space, flow, dt);
        ComputeFunctionG(gh, NULL, NULL, k + 1, j, i, U, space, flow, dt);
        ComputeNumericalDissipationDelta(delta, NULL, NULL, k, j, i, U, space, flow);
        CalculateGamma(gamma, g, gh, alpha);
        for (int row = 0; row < 5; ++row) {
            Phiz[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
        }
    }
    if (NULL != Phiy) {
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(NULL, lambda, NULL, NULL, alpha, NULL, k, j, i, U, space, flow);
        ComputeFunctionG(NULL, g, NULL, k, j, i, U, space, flow, dt);
        ComputeFunctionG(NULL, gh, NULL, k, j + 1, i, U, space, flow, dt);
        ComputeNumericalDissipationDelta(NULL, delta, NULL, k, j, i, U, space, flow);
        CalculateGamma(gamma, g, gh, alpha);
        for (int row = 0; row < 5; ++row) {
            Phiy[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
        }
    }
    if (NULL != Phix) {
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(NULL, NULL, lambda, NULL, NULL, alpha, k, j, i, U, space, flow);
        ComputeFunctionG(NULL, NULL, g, k, j, i, U, space, flow, dt);
        ComputeFunctionG(NULL, NULL, gh, k, j, i + 1, U, space, flow, dt);
        ComputeNumericalDissipationDelta(NULL, NULL, delta, k, j, i, U, space, flow);
        CalculateGamma(gamma, g, gh, alpha);
        for (int row = 0; row < 5; ++row) {
            Phix[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
        }
    }
    return 0;
}
static int ComputeFunctionG(
        Real gz[], Real gy[], Real gx[], 
        const int k, const int j, const int i, 
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
    if (NULL != gz) {
        const Real r = dt * space->ddz;
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(lambda, NULL, NULL, alpha, NULL, NULL, k, j, i, U, space, flow);
        ComputeNumericalDissipationDelta(delta, NULL, NULL, k, j, i, U, space, flow);
        CalculateSigma(sigma, lambda, delta, r);
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(lambdah, NULL, NULL, alphah, NULL, NULL, k - 1, j, i, U, space, flow);
        ComputeNumericalDissipationDelta(deltah, NULL, NULL, k - 1, j, i, U, space, flow);
        CalculateSigma(sigmah, lambdah, deltah, r);
        for (int row = 0; row < 5; ++row) {
            gz[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
        }
    }
    if (NULL != gy) {
        const Real r = dt * space->ddy;
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(NULL, lambda, NULL, NULL, alpha, NULL, k, j, i, U, space, flow);
        ComputeNumericalDissipationDelta(NULL, delta, NULL, k, j, i, U, space, flow);
        CalculateSigma(sigma, lambda, delta, r);
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(NULL, lambdah, NULL, NULL, alphah, NULL, k, j - 1, i, U, space, flow);
        ComputeNumericalDissipationDelta(NULL, deltah, NULL, k, j - 1, i, U, space, flow);
        CalculateSigma(sigmah, lambdah, deltah, r);
        for (int row = 0; row < 5; ++row) {
            gy[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
        }
    }
    if (NULL != gx) {
        const Real r = dt * space->ddx;
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(NULL, NULL, lambda, NULL, NULL, alpha, k, j, i, U, space, flow);
        ComputeNumericalDissipationDelta(NULL, NULL, delta, k, j, i, U, space, flow);
        CalculateSigma(sigma, lambda, delta, r);
        ComputeEigenvaluesAndDecompositionCoefficientAlpha(NULL, NULL, lambdah, NULL, NULL, alphah, k, j, i - 1, U, space, flow);
        ComputeNumericalDissipationDelta(NULL, NULL, deltah, k, j, i - 1, U, space, flow);
        CalculateSigma(sigmah, lambdah, deltah, r);
        for (int row = 0; row < 5; ++row) {
            gx[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
        }
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
static int ComputeEigenvaluesAndDecompositionCoefficientAlpha(
        Real lambdaz[], Real lambday[], Real lambdax[], 
        Real alphaz[], Real alphay[], Real alphax[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = ((k * space->jMax + j) * space->iMax + i) * 5;
    Real L[5][5] = {{0.0}}; /* store left eigenvectors */
    if (NULL != alphaz) {
        const int idxh = (((k + 1) * space->jMax + j) * space->iMax + i) * 5;
        const Real deltaU[5] = {
            U[idxh+0] - U[idx+0],
            U[idxh+1] - U[idx+1],
            U[idxh+2] - U[idx+2],
            U[idxh+3] - U[idx+3],
            U[idxh+4] - U[idx+4]};
        ComputeEigenvaluesAndEigenvectorSpaceL(lambdaz, NULL, NULL, L, NULL, NULL, k, j, i, U, space, flow);
        CalculateAlpha(alphaz, L, deltaU);
    }
    if (NULL != alphay) {
        const int idxh = ((k * space->jMax + j + 1) * space->iMax + i) * 5;
        const Real deltaU[5] = {
            U[idxh+0] - U[idx+0],
            U[idxh+1] - U[idx+1],
            U[idxh+2] - U[idx+2],
            U[idxh+3] - U[idx+3],
            U[idxh+4] - U[idx+4]};
        ComputeEigenvaluesAndEigenvectorSpaceL(NULL, lambday, NULL, NULL, L, NULL, k, j, i, U, space, flow);
        CalculateAlpha(alphay, L, deltaU);
    }
    if (NULL != alphax) {
        const int idxh = ((k * space->jMax + j) * space->iMax + i + 1) * 5;
        const Real deltaU[5] = {
            U[idxh+0] - U[idx+0],
            U[idxh+1] - U[idx+1],
            U[idxh+2] - U[idx+2],
            U[idxh+3] - U[idx+3],
            U[idxh+4] - U[idx+4]};
        ComputeEigenvaluesAndEigenvectorSpaceL(NULL, NULL, lambdax, NULL, NULL, L, k, j, i, U, space, flow);
        CalculateAlpha(alphax, L, deltaU);
    }
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
static int ComputeEigenvaluesAndEigenvectorSpaceL(
        Real lambdaz[], Real lambday[], Real lambdax[], 
        Real Lz[][5], Real Ly[][5], Real Lx[][5], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    if ((NULL != lambdaz) || (NULL != Lz)) {
        ComputeRoeAverage(Uo, NULL, NULL, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real c = Uo[5];
        const Real q = 0.5 * (u * u + v * v + w * w);
        const Real b = (flow->gamma - 1) / (2 * c * c);
        const Real d = (1 / (2 * c)); 
        if (NULL != lambdaz) {
            lambdaz[0] = w - c; lambdaz[1] = w; lambdaz[2] = w; lambdaz[3] = w; lambdaz[4] = w + c;
        }
        if (NULL != Lz) {
            Lz[0][0] = b * q + d * w;   Lz[0][1] = -b * u;             Lz[0][2] = -b * v;             Lz[0][3] = -b * w - d;     Lz[0][4] = b;
            Lz[1][0] = -2 * b * q * u;  Lz[1][1] = 2 * b * u * u + 1;  Lz[1][2] = 2 * b * v * u;      Lz[1][3] = 2 * b * w * u;  Lz[1][4] = -2 * b * u;
            Lz[2][0] = -2 * b * q * v;  Lz[2][1] = 2 * b * v * u;      Lz[2][2] = 2 * b * v * v + 1;  Lz[2][3] = 2 * b * w * v;  Lz[2][4] = -2 * b * v;
            Lz[3][0] = -2 * b * q + 1;  Lz[3][1] = 2 * b * u;          Lz[3][2] = 2 * b * v;          Lz[3][3] = 2 * b * w;      Lz[3][4] = -2 * b;
            Lz[4][0] = b * q - d * w;   Lz[4][1] = -b * u;             Lz[4][2] = -b * v;             Lz[4][3] = -b * w + d;     Lz[4][4] = b;
        }
    }
    if ((NULL != lambday) || (NULL != Ly)) {
        ComputeRoeAverage(NULL, Uo, NULL, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real c = Uo[5];
        const Real q = 0.5 * (u * u + v * v + w * w);
        const Real b = (flow->gamma - 1) / (2 * c * c);
        const Real d = (1 / (2 * c)); 
        if (NULL != lambday) {
            lambday[0] = v - c; lambday[1] = v; lambday[2] = v; lambday[3] = v; lambday[4] = v + c;
        }
        if (NULL != Ly) {
            Ly[0][0] = b * q + d * v;    Ly[0][1] = -b * u;             Ly[0][2] = -b * v - d;     Ly[0][3] = -b * w;             Ly[0][4] = b;
            Ly[1][0] = -2 * b * q * u;   Ly[1][1] = 2 * b * u * u + 1;  Ly[1][2] = 2 * b * v * u;  Ly[1][3] = 2 * b * w * u;      Ly[1][4] = -2 * b * u;
            Ly[2][0] = -2 * b * q + 1;   Ly[2][1] = 2 * b * u;          Ly[2][2] = 2 * b * v;      Ly[2][3] = 2 * b * w;          Ly[2][4] = -2 * b;
            Ly[3][0] = -2 * b * q * w;   Ly[3][1] = 2 * b * w * u;      Ly[3][2] = 2 * b * w * v;  Ly[3][3] = 2 * b * w * w + 1;  Ly[3][4] = -2 * b * w;
            Ly[4][0] = b * q - d * v;    Ly[4][1] = -b * u;             Ly[4][2] = -b * v + d;     Ly[4][3] = -b * w;             Ly[4][4] = b;
        }
    }
    if ((NULL != lambdax) || (NULL != Lx)) {
        ComputeRoeAverage(NULL, NULL, Uo, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real c = Uo[5];
        const Real q = 0.5 * (u * u + v * v + w * w);
        const Real b = (flow->gamma - 1) / (2 * c * c);
        const Real d = (1 / (2 * c)); 
        if (NULL != lambdax) {
            lambdax[0] = u - c; lambdax[1] = u; lambdax[2] = u; lambdax[3] = u; lambdax[4] = u + c;
        }
        if (NULL != Lx) {
            Lx[0][0] = b * q + d * u;   Lx[0][1] = -b * u - d;     Lx[0][2] = -b * v;             Lx[0][3] = -b * w;             Lx[0][4] = b;
            Lx[1][0] = -2 * b * q + 1;  Lx[1][1] = 2 * b * u;      Lx[1][2] = 2 * b * v;          Lx[1][3] = 2 * b * w;          Lx[1][4] = -2 * b;
            Lx[2][0] = -2 * b * q * v;  Lx[2][1] = 2 * b * v * u;  Lx[2][2] = 2 * b * v * v + 1;  Lx[2][3] = 2 * b * w * v;      Lx[2][4] = -2 * b * v;
            Lx[3][0] = -2 * b * q * w;  Lx[3][1] = 2 * b * w * u;  Lx[3][2] = 2 * b * w * v;      Lx[3][3] = 2 * b * w * w + 1;  Lx[3][4] = -2 * b * w;
            Lx[4][0] = b * q - d * u;   Lx[4][1] = -b * u + d;     Lx[4][2] = -b * v;             Lx[4][3] = -b * w;             Lx[4][4] = b;
        }
    }
    return 0;
}
static int ComputeEigenvectorSpaceR(
        Real Rz[][5], Real Ry[][5], Real Rx[][5], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    if (NULL != Rz) {
        ComputeRoeAverage(Uo, NULL, NULL, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real hT = Uo[4];
        const Real c = Uo[5];
        const Real q = 0.5 * (u * u + v * v + w * w);
        Rz[0][0] = 1;           Rz[0][1] = 0;  Rz[0][2] = 0;  Rz[0][3] = 1;          Rz[0][4] = 1;
        Rz[1][0] = u;           Rz[1][1] = 1;  Rz[1][2] = 0;  Rz[1][3] = 0;          Rz[1][4] = u;
        Rz[2][0] = v;           Rz[2][1] = 0;  Rz[2][2] = 1;  Rz[2][3] = 0;          Rz[2][4] = v;
        Rz[3][0] = w - c;       Rz[3][1] = 0;  Rz[3][2] = 0;  Rz[3][3] = w;          Rz[3][4] = w + c;
        Rz[4][0] = hT - w * c;  Rz[4][1] = u;  Rz[4][2] = v;  Rz[4][3] = w * w - q;  Rz[4][4] = hT + w * c;
    }
    if (NULL != Ry) {
        ComputeRoeAverage(NULL, Uo, NULL, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real hT = Uo[4];
        const Real c = Uo[5];
        const Real q = 0.5 * (u * u + v * v + w * w);
        Ry[0][0] = 1;           Ry[0][1] = 0;  Ry[0][2] = 1;          Ry[0][3] = 0;  Ry[0][4] = 1;
        Ry[1][0] = u;           Ry[1][1] = 1;  Ry[1][2] = 0;          Ry[1][3] = 0;  Ry[1][4] = u;
        Ry[2][0] = v - c;       Ry[2][1] = 0;  Ry[2][2] = v;          Ry[2][3] = 0;  Ry[2][4] = v + c;
        Ry[3][0] = w;           Ry[3][1] = 0;  Ry[3][2] = 0;          Ry[3][3] = 1;  Ry[3][4] = w;
        Ry[4][0] = hT - v * c;  Ry[4][1] = u;  Ry[4][2] = v * v - q;  Ry[4][3] = w;  Ry[4][4] = hT + v * c;
    }
    if (NULL != Rx) {
        ComputeRoeAverage(NULL, NULL, Uo, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real hT = Uo[4];
        const Real c = Uo[5];
        const Real q = 0.5 * (u * u + v * v + w * w);
        Rx[0][0] = 1;           Rx[0][1] = 1;          Rx[0][2] = 0;  Rx[0][3] = 0;  Rx[0][4] = 1;
        Rx[1][0] = u - c;       Rx[1][1] = u;          Rx[1][2] = 0;  Rx[1][3] = 0;  Rx[1][4] = u + c;
        Rx[2][0] = v;           Rx[2][1] = 0;          Rx[2][2] = 1;  Rx[2][3] = 0;  Rx[2][4] = v;
        Rx[3][0] = w;           Rx[3][1] = 0;          Rx[3][2] = 0;  Rx[3][3] = 1;  Rx[3][4] = w;
        Rx[4][0] = hT - u * c;  Rx[4][1] = u * u - q;  Rx[4][2] = v;  Rx[4][3] = w;  Rx[4][4] = hT + u * c;
    }
    return 0;
}
static int ComputeRoeAverage(
        Real Uoz[], Real Uoy[], Real Uox[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = ((k * space->jMax + j) * space->iMax + i) * 5;
    if (NULL != Uoz) {
        const int idxh = (((k + 1) * space->jMax + j) * space->iMax + i) * 5;
        CalculateRoeAverageUo(Uoz, idx, idxh, U, flow);
    }
    if (NULL != Uoy) {
        const int idxh = ((k * space->jMax + j + 1) * space->iMax + i) * 5;
        CalculateRoeAverageUo(Uoy, idx, idxh, U, flow);
    }
    if (NULL != Uox) {
        const int idxh = ((k * space->jMax + j) * space->iMax + i + 1) * 5;
        CalculateRoeAverageUo(Uox, idx, idxh, U, flow);
    }
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
    const Real hT = flow->gamma * U[idx+4] / rho - 0.5 * (u * u + v * v + w * w) * (flow->gamma - 1);
    const Real rho_h = U[idxh+0];
    const Real u_h = U[idxh+1] / rho_h;
    const Real v_h = U[idxh+2] / rho_h;
    const Real w_h = U[idxh+3] / rho_h;
    const Real hT_h = flow->gamma * U[idxh+4] / rho_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h) * (flow->gamma - 1);
    const Real D = sqrt(rho_h / rho);
    //Uo[0] = rho * (0.5 * (1 + D)) * (0.5 * (1 + D)); /* rho average, not required */
    Uo[1] = (u + D * u_h) / (1 + D); /* u average */
    Uo[2] = (v + D * v_h) / (1 + D); /* v average */
    Uo[3] = (w + D * w_h) / (1 + D); /* w average */
    Uo[4] = (hT + D * hT_h) / (1 + D); /* hT average */
    Uo[5] = sqrt((flow->gamma - 1) * (Uo[4] - 0.5 * (Uo[1] * Uo[1] + Uo[2] * Uo[2] + Uo[3] * Uo[3]))); /* the speed of sound */
    return 0;
}
static int ComputeNonViscousFlux(
        Real Fz[], Real Fy[], Real Fx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    const int idx = ((k * space->jMax + j) * space->iMax + i) * 5;
    const Real rho = U[idx+0];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = (flow->gamma - 1) * rho * (eT - 0.5 * (u * u + v * v + w * w));

    if (NULL != Fz) {
        Fz[0] = rho * w;
        Fz[1] = rho * w * u;
        Fz[2] = rho * w * v;
        Fz[3] = rho * w * w + p;
        Fz[4] = (rho * eT + p) * w;
    }
    if (NULL != Fy) {
        Fy[0] = rho * v;
        Fy[1] = rho * v * u;
        Fy[2] = rho * v * v + p;
        Fy[3] = rho * v * w;
        Fy[4] = (rho * eT + p) * v;
    }
    if (NULL != Fx) {
        Fx[0] = rho * u;
        Fx[1] = rho * u * u + p;
        Fx[2] = rho * u * v;
        Fx[3] = rho * u * w;
        Fx[4] = (rho * eT + p) * u;
    }
    return 0;
}
static int ComputeViscousFluxGradient(
        Real gradGz[], Real gradGy[], Real gradGx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
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
    Real hG[5] = {0.0}; /* viscous flux vector */
    Real Gh[5] = {0.0}; /* viscous flux vector */
    Real h = 0; /* reciprocal of differencing distance */
    if (NULL != gradGz) {
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
        const int idxl = (hl * space->jMax + j) * space->iMax + i;
        const int idxr = (hr * space->jMax + j) * space->iMax + i;
        if (10 <= space->nodeFlag[idxl]) {
            ++hl;
        }
        if (10 <= space->nodeFlag[idxr]) {
            --hr;
        }
        if (0 != (hr - hl)) { /* only do calculation when needed */
            h = space->ddz / (hr - hl);
            ComputeViscousFlux(hG, NULL, NULL, hl, j, i, U, space, flow);
            ComputeViscousFlux(Gh, NULL, NULL, hr, j, i, U, space, flow);
        } else {
            h = 0;
        }
        for (int row = 0; row < 5; ++row) {
            gradGz[row] = h * (Gh[row] - hG[row]);
        }
    }
    if (NULL != gradGy) {
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
        const int idxl = (k * space->jMax + hl) * space->iMax + i;
        const int idxr = (k * space->jMax + hr) * space->iMax + i;
        if (10 <= space->nodeFlag[idxl]) {
            ++hl;
        }
        if (10 <= space->nodeFlag[idxr]) {
            --hr;
        }
        if (0 != (hr - hl)) { /* only do calculation when needed */
            h = space->ddy / (hr - hl);
            ComputeViscousFlux(NULL, hG, NULL, k, hl, i, U, space, flow);
            ComputeViscousFlux(NULL, Gh, NULL, k, hr, i, U, space, flow);
        } else {
            h = 0;
        }
        for (int row = 0; row < 5; ++row) {
            gradGy[row] = h * (Gh[row] - hG[row]);
        }
    }
    if (NULL != gradGx) {
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
        const int idxl = (k * space->jMax + j) * space->iMax + hl;
        const int idxr = (k * space->jMax + j) * space->iMax + hr;
        if (10 <= space->nodeFlag[idxl]) {
            ++hl;
        }
        if (10 <= space->nodeFlag[idxr]) {
            --hr;
        }
        if (0 != (hr - hl)) { /* only do calculation when needed */
            h = space->ddx / (hr - hl);
            ComputeViscousFlux(NULL, NULL, hG, k, j, hl, U, space, flow);
            ComputeViscousFlux(NULL, NULL, Gh, k, j, hr, U, space, flow);
        } else {
            h = 0;
        }
        for (int row = 0; row < 5; ++row) {
            gradGx[row] = h * (Gh[row] - hG[row]);
        }
    }
    return 0;
}
static int ComputeViscousFlux(
        Real Gz[], Real Gy[], Real Gx[], 
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow)
{
    Real rho = 0; 
    Real rho_h = 0; 
    Real u = 0;
    Real u_h = 0;
    Real v = 0;
    Real v_h = 0;
    Real w = 0;
    Real w_h = 0;
    Real eT = 0;
    Real eT_h = 0;
    Real T = 0;
    Real T_h = 0;
    const int idx = ((k * space->jMax + j) * space->iMax + i) * 5;
    const int idxW = ((k * space->jMax + j) * space->iMax + i - 1) * 5;
    const int idxE = ((k * space->jMax + j) * space->iMax + i + 1) * 5;
    const int idxS = ((k * space->jMax + j - 1) * space->iMax + i) * 5;
    const int idxN = ((k * space->jMax + j + 1) * space->iMax + i) * 5;
    const int idxF = (((k - 1) * space->jMax + j) * space->iMax + i) * 5;
    const int idxB = (((k + 1) * space->jMax + j) * space->iMax + i) * 5;

    /* calculate derivatives in z direction */
    rho_h = U[idxB+0];
    u_h = U[idxB+1] / rho_h;
    v_h = U[idxB+2] / rho_h;
    w_h = U[idxB+3] / rho_h;
    eT_h = U[idxB+4] / rho_h;
    T_h = (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h)) / flow->cv;

    rho = U[idxF+0];
    u = U[idxF+1] / rho;
    v = U[idxF+2] / rho;
    w = U[idxF+3] / rho;
    eT = U[idxF+4] / rho;
    T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

    const Real du_dz = (u_h - u) * (0.5 * space->ddz);
    const Real dv_dz = (v_h - v) * (0.5 * space->ddz);
    const Real dw_dz = (w_h - w) * (0.5 * space->ddz);
    const Real dT_dz = (T_h - T) * (0.5 * space->ddz);

    /* calculate derivatives in y direction */
    rho_h = U[idxN+0];
    u_h = U[idxN+1] / rho_h;
    v_h = U[idxN+2] / rho_h;
    w_h = U[idxN+3] / rho_h;
    eT_h = U[idxN+4] / rho_h;
    T_h = (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h)) / flow->cv;

    rho = U[idxS+0];
    u = U[idxS+1] / rho;
    v = U[idxS+2] / rho;
    w = U[idxS+3] / rho;
    eT = U[idxS+4] / rho;
    T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

    const Real du_dy = (u_h - u) * (0.5 * space->ddy);
    const Real dv_dy = (v_h - v) * (0.5 * space->ddy);
    const Real dw_dy = (w_h - w) * (0.5 * space->ddy);
    const Real dT_dy = (T_h - T) * (0.5 * space->ddy);

    /* calculate derivatives in x direction */
    rho_h = U[idxE+0];
    u_h = U[idxE+1] / rho_h;
    v_h = U[idxE+2] / rho_h;
    w_h = U[idxE+3] / rho_h;
    eT_h = U[idxE+4] / rho_h;
    T_h = (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h)) / flow->cv;

    rho = U[idxW+0];
    u = U[idxW+1] / rho;
    v = U[idxW+2] / rho;
    w = U[idxW+3] / rho;
    eT = U[idxW+4] / rho;
    T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

    const Real du_dx = (u_h - u) * (0.5 * space->ddx);
    const Real dv_dx = (v_h - v) * (0.5 * space->ddx);
    const Real dw_dx = (w_h - w) * (0.5 * space->ddx);
    const Real dT_dx = (T_h - T) * (0.5 * space->ddx);

    /* regain the primitive variables in current point */
    rho = U[idx+0];
    u = U[idx+1] / rho;
    v = U[idx+2] / rho;
    w = U[idx+3] / rho;
    eT = U[idx+4] / rho;
    T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

    /* Calculate dynamic viscosity and heat conductivity */
    const Real mu = flow->refMu * 1.45e-6 * (pow(T * flow->refTemperature, 1.5) / (T * flow->refTemperature + 110));
    const Real heatK = flow->gamma * flow->cv * mu / flow->refPr;
    const Real divV = du_dx + dv_dy + dw_dz;

    if (NULL != Gz) {
        Gz[0] = 0;
        Gz[1] = mu * (dw_dx + du_dz);
        Gz[2] = mu * (dw_dy + dv_dz);
        Gz[3] = mu * (2 * dw_dz - (2/3) * divV);
        Gz[4] = heatK * dT_dz + 
            u * Gz[1] + v * Gz[2] + w * Gz[3];
    }
    if (NULL != Gy) {
        Gy[0] = 0;
        Gy[1] = mu * (dv_dx + du_dy);
        Gy[2] = mu * (2 * dv_dy - (2/3) * divV);
        Gy[3] = mu * (dv_dz + dw_dy);
        Gy[4] = heatK * dT_dy + 
            u * Gy[1] + v * Gy[2] + w * Gy[3];
    }
    if (NULL != Gx) {
        Gx[0] = 0;
        Gx[1] = mu * (2 * du_dx - (2/3) * divV);
        Gx[2] = mu * (du_dy + dv_dx);
        Gx[3] = mu * (du_dz + dw_dx);
        Gx[4] = heatK * dT_dx + 
            u * Gx[1] + v * Gx[2] + w * Gx[3];
    }
    return 0;
}
static int ComputeNumericalDissipationDelta(
        Real deltaz[], Real deltay[], Real deltax[], 
        const int k, const int j, const int i,
        const Real *U, const Space *space, const Flow *flow)
{
    Real Uo[6] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    const Real delta0 = 0.125; /* in [0.05, 0.25], 0.125 is recommended */
    if (NULL != deltaz) {
        ComputeRoeAverage(Uo, NULL, NULL, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real c = Uo[5];
        for (int row = 0; row < 5; ++row) {
            deltaz[row] = delta0 * (fabs(u) + fabs(v) + fabs(w) + c); 
        }
    }
    if (NULL != deltay) {
        ComputeRoeAverage(NULL, Uo, NULL, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real c = Uo[5];
        for (int row = 0; row < 5; ++row) {
            deltay[row] = delta0 * (fabs(u) + fabs(v) + fabs(w) + c); 
        }
    }
    if (NULL != deltax) {
        ComputeRoeAverage(NULL, NULL, Uo, k, j, i, U, space, flow);
        const Real u = Uo[1];
        const Real v = Uo[2];
        const Real w = Uo[3];
        const Real c = Uo[5];
        for (int row = 0; row < 5; ++row) {
            deltax[row] = delta0 * (fabs(u) + fabs(v) + fabs(w) + c); 
        }
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
/* a good practice: end file with a newline */

