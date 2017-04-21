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
 * Required Header Files
 ****************************************************************************/
#include "solid_dynamics.h"
#include <stdlib.h> /* mathematical functions on integers */
#include <math.h> /* common mathematical functions */
#include <limits.h> /* sizes of integral types */
#include <float.h> /* size of floating point values */
#include <string.h> /* manipulating strings */
#include "immersed_boundary.h"
#include "computational_geometry.h"
#include "linear_system.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void Kinematics(const Real, const Real, Space *);
static void CollisionDynamics(Space *);
static void CollisionState(const int, const int, const int, const int, const int,
        const int [restrict][DIMS], const Node *const, const Partition *, Geometry *);
static void AddCollideObject(const int [restrict], const int, Geometry *);
static void ApplyMotion(const Real, Space *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int SolidDynamics(const Real now, const Real dt, Space *space, const Model *model)
{
    /*
     * Compute exerted surface forces.
     */
    SurfaceForceIntegration(space, model);
    /*
     * Kinematics
     */
    Kinematics(now, dt, space);
    /*
     * Collision
     */
    if (1 != model->fsi) {
        CollisionDynamics(space);
    }
    /*
     * Update spatial position
     */
    ApplyMotion(dt, space);
    /*
     * Recompute the geometry domain.
     */
    ComputeGeometryDomain(space, model);
    /*
     * Apply boundary condition
     */
    ImmersedBoundaryTreatment(TO, space, model);
    return 0;
}
void SurfaceForceIntegration(Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    const Node *const node = space->node;
    const IntVec nMin = {part->ns[PIN][X][MIN], part->ns[PIN][Y][MIN], part->ns[PIN][Z][MIN]};
    const IntVec nMax = {part->ns[PIN][X][MAX], part->ns[PIN][Y][MAX], part->ns[PIN][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    int idx = 0; /* linear array index math variable */
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    int lidN = 0; /* count total number of interfacial nodes */
    int gstN = 0; /* count total number of ghost nodes */
    const Real zero = 0.0;
    const Real percent = FLT_EPSILON * FLT_EPSILON;
    RealVec pG = {zero}; /* ghost point */
    RealVec pO = {zero}; /* boundary point */
    RealVec pI = {zero}; /* image point */
    RealVec N = {zero}; /* normal */
    Real Uo[DIMUo] = {zero};
    RealVec V = {zero}; /* velocity vector */
    RealVec r = {zero}; /* position vector */
    RealVec Fp = {zero}; /* pressure force */
    RealVec Fv = {zero}; /* viscous force */
    RealVec Fs = {zero}; /* surface force */
    RealVec Tt = {zero}; /* torque */
    RealVec fvar = {zero}; /* force offset, mean, variance */
    Real Vn = zero; /* velocity projection */
    Real mu = zero; /* viscosity */
    Real ds = zero; /* infinitesimal area for integration */
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        if (0 < poly->state) { /* surface force negligible */
            continue;
        }
        /* reset some non accumulative information to zero */
        memset(poly->Fp, 0, DIMS * sizeof(*poly->Fp));
        memset(poly->Fv, 0, DIMS * sizeof(*poly->Fv));
        memset(poly->Tt, 0, DIMS * sizeof(*poly->Tt));
        memset(fvar, 0, DIMS * sizeof(*fvar));
        lidN = 0;
        gstN = 0;
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ValidNodeSpace(NodeSpace(poly->box[s][MIN], sMin[s], dd[s], ng), nMin[s], nMax[s]);
            box[s][MAX] = ValidNodeSpace(NodeSpace(poly->box[s][MAX], sMin[s], dd[s], ng), nMin[s], nMax[s]) + 1;
        }
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if ((2 == node[idx].lid) && (n + 1 == node[idx].gid)) {
                        ++lidN; /* an interfacial node of current geometry */
                    }
                    if ((2 != node[idx].gst) || (n + 1 != node[idx].gid)) {
                        continue;
                    }
                    ++gstN; /* a ghost node of current geometry */
                    /* surface force exerted by fluid (pressure + shear force) */
                    pG[X] = PointSpace(i, sMin[X], d[X], ng);
                    pG[Y] = PointSpace(j, sMin[Y], d[Y], ng);
                    pG[Z] = PointSpace(k, sMin[Z], d[Z], ng);
                    ComputeGeometricData(node[idx].fid, poly, pG, pO, pI, N);
                    r[X] = pO[X] - poly->O[X];
                    r[Y] = pO[Y] - poly->O[Y];
                    r[Z] = pO[Z] - poly->O[Z];
                    PrimitiveByConservative(model->gamma, model->gasR, node[idx].U[TO], Uo);
                    Fp[X] = Uo[4] * N[X];
                    Fp[Y] = Uo[4] * N[Y];
                    Fp[Z] = Uo[4] * N[Z];
                    if (1 == gstN) {
                        fvar[0] = Uo[4];
                    }
                    fvar[1] = fvar[1] + Uo[4] - fvar[0];
                    fvar[2] = fvar[2] + (Uo[4] - fvar[0]) * (Uo[4] - fvar[0]);
                    if ((zero < model->refMu) && (zero < poly->cf)) {
                        mu = model->refMu * Viscosity(Uo[5] * model->refT);
                        Cross(poly->W[TO], r, V);
                        V[X] = Uo[1] - (poly->V[TO][X] + V[X]);
                        V[Y] = Uo[2] - (poly->V[TO][Y] + V[Y]);
                        V[Z] = Uo[3] - (poly->V[TO][Z] + V[Z]);
                        Vn = Dot(V, N);
                        Fv[X] = mu * (V[X] - Vn * N[X]) / Dist(pG, pO);
                        Fv[Y] = mu * (V[Y] - Vn * N[Y]) / Dist(pG, pO);
                        Fv[Z] = mu * (V[Z] - Vn * N[Z]) / Dist(pG, pO);
                    } else {
                        memset(Fv, 0, DIMS * sizeof(*Fv));
                    }
                    Fs[X] = Fp[X] + Fv[X];
                    Fs[Y] = Fp[Y] + Fv[Y];
                    Fs[Z] = Fp[Z] + Fv[Z];
                    Cross(r, Fs, Tt);
                    /* integration sum */
                    for (int s = 0; s < DIMS; ++s) {
                        poly->Fp[s] = poly->Fp[s] + Fp[s];
                        poly->Fv[s] = poly->Fv[s] + Fv[s];
                        poly->Tt[s] = poly->Tt[s] + Tt[s];
                    }
                }
            }
        }
        /* calibrate the sum of discrete forces into integration */
        if ((0 == lidN) || (0 == gstN)) { /* no surface force exerted */
            continue;
        }
        ds = poly->area / lidN;
        fvar[2] = (fvar[2] - fvar[1] * fvar[1] / gstN) / gstN; /* variance */
        fvar[1] = fvar[1] / gstN + fvar[0]; /* mean */
        if (percent * fvar[1] * fvar[1] > fvar[2]) { /* recover equilibrium state and ignore integration error */
            ds = zero;
        }
        for (int s = 0; s < DIMS; ++s) {
            poly->Fp[s] = -poly->Fp[s] * ds;
            poly->Fv[s] = -poly->Fv[s] * ds;
            poly->Tt[s] = -poly->Tt[s] * ds;
        }
    }
    return;
}
static void Kinematics(const Real now, const Real dt, Space *space)
{
    Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    Real A[DIMS][DIMS] = {{0.0}};
    Real B[DIMS][1] = {{0.0}};
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        if (1 == poly->state) { /* stationary object */
            continue;
        }
        if (now > poly->to) { /* end power supply */
            memset(poly->at[TN], 0, DIMS * sizeof(*poly->at[TN]));
            memset(poly->ar[TN], 0, DIMS * sizeof(*poly->ar[TN]));
            poly->to = FLT_MAX; /* avoid repeating */
        }
        /* translation and rotational acceleration */
        for (int s = 0; s < DIMS; ++s) {
            for (int m = 0; m < DIMS; ++m) {
                A[s][m] = poly->I[s][m];
            }
            B[s][0] = poly->Tt[s] / poly->rho;
        }
        MatrixLinearSystemSolver(DIMS, A, 1, B, B);
        for (int s = 0; s < DIMS; ++s) {
            /* acceleration from surface force and body force */
            poly->at[TO][s] = (poly->Fp[s] + poly->Fv[s]) / (poly->rho * poly->volume) + poly->at[TN][s] + poly->g[s];
            poly->ar[TO][s] = B[s][0] + poly->ar[TN][s];
        }
        /* velocity integration */
        for (int s = 0; s < DIMS; ++s) {
            /* averaged velocity during time level n and n+1 */
            poly->V[TN][s] = 0.5 * (poly->V[TO][s] + poly->V[TO][s] + poly->at[TO][s] * dt);
            poly->W[TN][s] = 0.5 * (poly->W[TO][s] + poly->W[TO][s] + poly->ar[TO][s] * dt);
            /* velocity at the next time level */
            poly->V[TO][s] = poly->V[TO][s] + poly->at[TO][s] * dt;
            poly->W[TO][s] = poly->W[TO][s] + poly->ar[TO][s] * dt;
        }
    }
    return;
}
static void CollisionDynamics(Space *space)
{
    const Partition *restrict part = &(space->part);
    Geometry *geo = &(space->geo);
    Collision *col = NULL;
    Polyhedron *polp = NULL;
    Polyhedron *poln = NULL;
    const Node *const node = space->node;
    const IntVec nMin = {part->ns[PIN][X][MIN], part->ns[PIN][Y][MIN], part->ns[PIN][Z][MIN]};
    const IntVec nMax = {part->ns[PIN][X][MAX], part->ns[PIN][Y][MAX], part->ns[PIN][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    int idx = 0; /* linear array index math variable */
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    const Real zero = 0.0;
    const Real one = 1.0;
    const int coltag = INT_MAX / 2; /* colliding polyhedron marker */
    const Real crList[5] = {0.0, 0.25, 0.5, 0.75, 1.0}; /* coefficient of restitution */
    RealVec Vo = {zero}; /* original translational velocity */
    RealVec Wo = {zero}; /* original rotational velocity */
    RealVec V = {zero}; /* relative translational velocity */
    RealVec W = {zero}; /* relative rotational velocity */
    RealVec N = {zero}; /* line of impact */
    Real Vn = zero; /* translational velocity projection on line of impact */
    Real cr = zero; /* coefficient of restitution */
    Real cf = zero; /* coefficient of sliding friction */
    Real mp = zero; /* mass */
    Real mn = zero; /* mass */
    Real meff = zero; /* effective mass */
    for (int p = 0; p < geo->totN; ++p) {
        polp = geo->poly + p;
        if (1 == polp->state) { /* stationary object */
            continue;
        }
        geo->colN = 0; /* reset */
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ValidNodeSpace(NodeSpace(polp->box[s][MIN], sMin[s], dd[s], ng), nMin[s], nMax[s]);
            box[s][MAX] = ValidNodeSpace(NodeSpace(polp->box[s][MAX], sMin[s], dd[s], ng), nMin[s], nMax[s]) + 1;
        }
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if ((1 != node[idx].lid) || (p + 1 != node[idx].gid)) {
                        continue;
                    }
                    CollisionState(k, j, i, p + 1, part->pathSep[1], part->path, node, part, geo);
                }
            }
        }
        /* skip none contacting polyhedron */
        if (0 == geo->colN) {
            continue;
        }
        /* backup original velocity */
        memcpy(Vo, polp->V[TO], DIMS * sizeof(*polp->V[TO]));
        memcpy(Wo, polp->W[TO], DIMS * sizeof(*polp->W[TO]));
        /* initialize post-collision velocity */
        memset(polp->V[TO], 0, DIMS * sizeof(*polp->V[TO]));
        memset(polp->W[TO], 0, DIMS * sizeof(*polp->W[TO]));
        /* pairwise collision */
        mp = polp->rho * polp->volume;
        for (int n = 0; n < geo->colN; ++n) {
            col = geo->col + n;
            if (0 == abs(col->N[X]) + abs(col->N[Y]) + abs(col->N[Z])) {
                if ((geo->colN - 1 == n) && (coltag > polp->state)) {
                    /* recover contacting but none colliding polyhedron */
                    memcpy(polp->V[TO], Vo, DIMS * sizeof(*polp->V[TO]));
                    memcpy(polp->W[TO], Wo, DIMS * sizeof(*polp->W[TO]));
                }
                continue;
            }
            N[X] = col->N[X];
            N[Y] = col->N[Y];
            N[Z] = col->N[Z];
            Normalize(DIMS, Norm(N), N);
            poln = geo->poly + col->gid - 1;
            mn = poln->rho * poln->volume;
            meff = mn / (mp + mn);
            /* relative speed */
            for (int s = 0; s < DIMS; ++s) {
                V[s] = polp->V[TN][s] - poln->V[TN][s];
                W[s] = polp->W[TN][s] - poln->W[TN][s];
            }
            cr = 0.5 * (crList[polp->mid] + crList[poln->mid]);
            cf = 0.5 * (polp->cf + poln->cf);
            Vn = Dot(V, N);
            if (zero >= Vn) {
                if ((geo->colN - 1 == n) && (coltag > polp->state)) {
                    /* recover contacting but none colliding polyhedron */
                    memcpy(polp->V[TO], Vo, DIMS * sizeof(*polp->V[TO]));
                    memcpy(polp->W[TO], Wo, DIMS * sizeof(*polp->W[TO]));
                }
                continue;
            }
            /* mark colliding polyhedron */
            if (coltag > polp->state) {
                polp->state = polp->state + coltag;
            }
            /* vector summation of the collision velocities in the global frame */
            for (int s = 0; s < DIMS; ++s) {
                polp->V[TO][s] = polp->V[TO][s] + polp->V[TN][s] - meff * (one + cr) * Vn * N[s] - cf * (V[s] - Vn * N[s]);
                polp->W[TO][s] = polp->W[TO][s] + polp->W[TN][s] - meff * W[s];
            }
        }
    }
    /* update post-collision velocity for collided polyhedron */
    for (int p = 0; p < geo->totN; ++p) {
        polp = geo->poly + p;
        if (coltag > polp->state) { /* none collided polyhedron */
            continue;
        }
        polp->state = polp->state - coltag; /* recover state */
        memcpy(polp->V[TN], polp->V[TO], DIMS * sizeof(*polp->V[TO]));
        memcpy(polp->W[TN], polp->W[TO], DIMS * sizeof(*polp->W[TO]));
    }
    return;
}
static void CollisionState(const int k, const int j, const int i, const int gid, const int end,
        const int path[restrict][DIMS], const Node *const node, const Partition *part, Geometry *geo)
{
    /*
     * Search around the specified node to find colliding objects.
     */
    int idx = 0; /* linear array index math variable */
    for (int n = 0; n < end; ++n) {
        idx = IndexNode(k + path[n][Z], j + path[n][Y], i + path[n][X], part->n[Y], part->n[X]);
        if (NONE == node[idx].gid) { /* an exterior node is not valid */
            continue;
        }
        if (0 == node[idx].gid) { /* a fluid node is not valid */
            continue;
        }
        if (gid != node[idx].gid) { /* a heterogeneous node on the path */
            AddCollideObject(path[n], node[idx].gid, geo);
        }
    }
    return;
}
static void AddCollideObject(const int N[restrict], const int gid, Geometry *geo)
{
    Collision *col = NULL;
    /* search the object list, if already exist, adjust the line of impact */
    for (int n = 0; n < geo->colN; ++n) {
        col = geo->col + n;
        if (gid == col->gid) {
            col->N[X] = col->N[X] + N[X];
            col->N[Y] = col->N[Y] + N[Y];
            col->N[Z] = col->N[Z] + N[Z];
            return;
        }
    }
    /* otherwise, add to the collision list */
    col = geo->col + geo->colN;
    col->gid = gid;
    col->N[X] = N[X];
    col->N[Y] = N[Y];
    col->N[Z] = N[Z];
    ++(geo->colN);
    return;
}
static void ApplyMotion(const Real dt, Space *space)
{
    Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    RealVec offset = {0.0}; /* translation */
    RealVec angle = {0.0}; /* rotation */
    const RealVec scale = {1.0, 1.0, 1.0}; /* scale */
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        if (1 == poly->state) { /* stationary object */
            continue;
        }
        /* position integration */
        for (int s = 0; s < DIMS; ++s) {
            offset[s] = poly->V[TN][s] * dt;
            angle[s] = poly->W[TN][s] * dt;
        }
        /* transform geometry */
        if (0 == poly->faceN) { /* analytical sphere */
            poly->O[X] = poly->O[X] + offset[X];
            poly->O[Y] = poly->O[Y] + offset[Y];
            poly->O[Z] = poly->O[Z] + offset[Z];
            /* bounding box */
            for (int s = 0; s < DIMS; ++s) {
                poly->box[s][MIN] = poly->O[s] - poly->r;
                poly->box[s][MAX] = poly->O[s] + poly->r;
            }
        } else { /* triangulated polyhedron */
            Transformation(poly->O, scale, angle, offset, poly);
        }
    }
    return;
}
/* a good practice: end file with a newline */

