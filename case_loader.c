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
#include "case_loader.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void ReadCaseSettingData(Time *, Space *, Model *);
static void ReadGeometrySettingData(Geometry *const);
static void ReadBoundaryData(FILE *, Space *, const int);
static void ReadConsecutiveData(FILE *, const int, const char *, Real *, char [][VARSTR]);
static void WriteBoundaryData(FILE *, const Space *, const int);
static void WriteInitializerData(FILE *, const Space *, const int);
static void WriteVerifyData(const Time *, const Space *, const Model *);
static void CheckCaseSettingData(const Time *, const Space *, const Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void LoadCaseData(Time *time, Space *space, Model *model)
{
    ReadCaseSettingData(time, space, model);
    ReadGeometrySettingData(&(space->geo));
    WriteVerifyData(time, space, model);
    CheckCaseSettingData(time, space, model);
    return;
}
static void ReadCaseSettingData(Time *time, Space *space, Model *model)
{
    Partition *const part = &(space->part);
    part->typeBC = AssignStorage(NBC * sizeof(*part->typeBC));
    part->N = AssignStorage(NBC * sizeof(*part->N));
    part->varBC = AssignStorage(NBC * sizeof(*part->varBC));
    part->typeIC = AssignStorage(NIC * sizeof(*part->typeIC));
    part->posIC = AssignStorage(NIC * sizeof(*part->posIC));
    part->varIC = AssignStorage(NIC * sizeof(*part->varIC));
    const char *fname = "artracfd.case";
    FILE *fp = Fopen(fname, "r");
    String str = {'\0'}; /* store the current read line */
    int nentry = 0; /* entry count */
    const char *fmtI = ParseFormat("%lg");
    const char *fmtJ = ParseFormat("%lg, %lg, %lg");
    while (NULL != fgets(str, sizeof str, fp)) {
        ParseCommand(str);
        if (0 == strncmp(str, "space begin", sizeof str)) {
            ++nentry;
            Sread(fp, 3, fmtJ, &(part->domain[X][MIN]), &(part->domain[Y][MIN]),
                    &(part->domain[Z][MIN]));
            Sread(fp, 3, fmtJ, &(part->domain[X][MAX]), &(part->domain[Y][MAX]),
                    &(part->domain[Z][MAX]));
            Sread(fp, 3, "%d, %d, %d", &(part->m[X]), &(part->m[Y]), &(part->m[Z]));
            continue;
        }
        if (0 == strncmp(str, "time begin", sizeof str)) {
            ++nentry;
            Sread(fp, 1, "%d", &(time->restart));
            Sread(fp, 1, fmtI, &(time->end));
            Sread(fp, 1, fmtI, &(time->numCFL));
            Sread(fp, 1, "%d", &(time->stepN));
            Sread(fp, 1, "%d", &(time->dataW[PROSD]));
            Sread(fp, 1, "%d", &(time->dataStreamer));
            continue;
        }
        if (0 == strncmp(str, "numerical begin", sizeof str)) {
            ++nentry;
            Sread(fp, 1, "%d", &(model->tScheme));
            Sread(fp, 1, "%d", &(model->sScheme));
            Sread(fp, 1, "%d", &(model->multidim));
            Sread(fp, 1, "%d", &(model->jacobMean));
            Sread(fp, 1, "%d", &(model->fluxSplit));
            Sread(fp, 1, "%d", &(model->psi));
            Sread(fp, 1, "%d", &(model->ibmLayer));
            continue;
        }
        if (0 == strncmp(str, "material begin", sizeof str)) {
            ++nentry;
            Sread(fp, 1, "%d", &(model->mid));
            Sread(fp, 1, fmtI, &(model->refMu));
            Sread(fp, 1, "%d", &(model->gState));
            Sread(fp, 3, fmtJ, &(model->g[X]), &(model->g[Y]), &(model->g[Z]));
            continue;
        }
        if (0 == strncmp(str, "reference begin", sizeof str)) {
            ++nentry;
            Sread(fp, 1, fmtI, &(model->refL));
            Sread(fp, 1, fmtI, &(model->refRho));
            Sread(fp, 1, fmtI, &(model->refV));
            Sread(fp, 1, fmtI, &(model->refT));
            continue;
        }
        if (0 == strncmp(str, "initialization begin", sizeof str)) {
            ++nentry;
            part->nIC = 0; /* enforce global initialization first */
            part->typeIC[part->nIC] = ICGLOBAL;
            ReadConsecutiveData(fp, VARIC, "%s", NULL, part->varIC[part->nIC]);
            ++part->nIC;
            continue;
        }
        if (0 == strncmp(str, "west boundary begin", sizeof str)) {
            ++nentry;
            ReadBoundaryData(fp, space, PWB);
            continue;
        }
        if (0 == strncmp(str, "east boundary begin", sizeof str)) {
            ++nentry;
            ReadBoundaryData(fp, space, PEB);
            continue;
        }
        if (0 == strncmp(str, "south boundary begin", sizeof str)) {
            ++nentry;
            ReadBoundaryData(fp, space, PSB);
            continue;
        }
        if (0 == strncmp(str, "north boundary begin", sizeof str)) {
            ++nentry;
            ReadBoundaryData(fp, space, PNB);
            continue;
        }
        if (0 == strncmp(str, "front boundary begin", sizeof str)) {
            ++nentry;
            ReadBoundaryData(fp, space, PFB);
            continue;
        }
        if (0 == strncmp(str, "back boundary begin", sizeof str)) {
            ++nentry;
            ReadBoundaryData(fp, space, PBB);
            continue;
        }
        if (0 == strncmp(str, "plane initialization begin", sizeof str)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->nIC] = ICPLANE;
            Sread(fp, 3, fmtJ, part->posIC[part->nIC] + 0,
                    part->posIC[part->nIC] + 1, part->posIC[part->nIC] + 2);
            Sread(fp, 3, fmtJ, part->posIC[part->nIC] + 3,
                    part->posIC[part->nIC] + 4, part->posIC[part->nIC] + 5);
            ReadConsecutiveData(fp, VARIC, "%s", NULL, part->varIC[part->nIC]);
            ++part->nIC;
            continue;
        }
        if (0 == strncmp(str, "sphere initialization begin", sizeof str)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->nIC] = ICSPHERE;
            Sread(fp, 3, fmtJ, part->posIC[part->nIC] + 0,
                    part->posIC[part->nIC] + 1, part->posIC[part->nIC] + 2);
            Sread(fp, 1, fmtI, part->posIC[part->nIC] + 6);
            ReadConsecutiveData(fp, VARIC, "%s", NULL, part->varIC[part->nIC]);
            ++part->nIC;
            continue;
        }
        if (0 == strncmp(str, "box initialization begin", sizeof str)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->nIC] = ICBOX;
            Sread(fp, 3, fmtJ, part->posIC[part->nIC] + 0,
                    part->posIC[part->nIC] + 1, part->posIC[part->nIC] + 2);
            Sread(fp, 3, fmtJ, part->posIC[part->nIC] + 3,
                    part->posIC[part->nIC] + 4, part->posIC[part->nIC] + 5);
            ReadConsecutiveData(fp, VARIC, "%s", NULL, part->varIC[part->nIC]);
            ++part->nIC;
            continue;
        }
        if (0 == strncmp(str, "cylinder initialization begin", sizeof str)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->nIC] = ICCYLINDER;
            Sread(fp, 3, fmtJ, part->posIC[part->nIC] + 0,
                    part->posIC[part->nIC] + 1, part->posIC[part->nIC] + 2);
            Sread(fp, 3, fmtJ, part->posIC[part->nIC] + 3,
                    part->posIC[part->nIC] + 4, part->posIC[part->nIC] + 5);
            Sread(fp, 1, fmtI, part->posIC[part->nIC] + 6);
            ReadConsecutiveData(fp, VARIC, "%s", NULL, part->varIC[part->nIC]);
            ++part->nIC;
            continue;
        }
        if (0 == strncmp(str, "probe count begin", sizeof str)) {
            /* optional entry do not increase entry count */
            Sread(fp, 1, "%d", &(time->dataN[PROPT]));
            Sread(fp, 1, "%d", &(time->dataN[PROLN]));
            Sread(fp, 1, "%d", &(time->dataN[PROCV]));
            Sread(fp, 1, "%d", &(time->dataN[PROFC]));
            if (0 < time->dataN[PROPT]) {
                time->pp = AssignStorage(time->dataN[PROPT] * sizeof(*time->pp));
            }
            if (0 < time->dataN[PROLN]) {
                time->lp = AssignStorage(time->dataN[PROLN] * sizeof(*time->lp));
            }
            continue;
        }
        if (0 == strncmp(str, "probe control begin", sizeof str)) {
            /* optional entry do not increase entry count */
            Sread(fp, 1, "%d", &(time->dataW[PROPT]));
            Sread(fp, 1, "%d", &(time->dataW[PROLN]));
            Sread(fp, 1, "%d", &(time->dataW[PROCV]));
            Sread(fp, 1, "%d", &(time->dataW[PROFC]));
            continue;
        }
        if (0 == strncmp(str, "point probe begin", sizeof str)) {
            /* optional entry do not increase entry count */
            for (int n = 0; n < time->dataN[PROPT]; ++n) {
                Sread(fp, 3, fmtJ, time->pp[n] + 0,
                        time->pp[n] + 1, time->pp[n] + 2);
            }
            continue;
        }
        if (0 == strncmp(str, "line probe begin", sizeof str)) {
            /* optional entry do not increase entry count */
            for (int n = 0; n < time->dataN[PROLN]; ++n) {
                Sread(fp, 3, fmtJ, time->lp[n] + 0,
                        time->lp[n] + 1, time->lp[n] + 2);
                Sread(fp, 3, fmtJ, time->lp[n] + 3,
                        time->lp[n] + 4, time->lp[n] + 5);
                Sread(fp, 1, fmtI, time->lp[n] + 6);
            }
            continue;
        }
    }
    fclose(fp);
    if (12 != nentry) {
        ShowError("missing or repeated sections: %s, entry: %d", fname, nentry);
    }
    return;
}
static void ReadGeometrySettingData(Geometry *const geo)
{
    const char *fname = "artracfd.geo";
    FILE *fp = Fopen(fname, "r");
    String str = {'\0'}; /* store the current read line */
    int nentry = 0; /* entry count */
    while (NULL != fgets(str, sizeof str, fp)) {
        ParseCommand(str);
        if (0 == strncmp(str, "count begin", sizeof str)) {
            ++nentry;
            Sread(fp, 1, "%d", &(geo->sphN));
            Sread(fp, 1, "%d", &(geo->stlN));
            break;
        }
    }
    fclose(fp);
    if (1 != nentry) {
        ShowError("missing or repeated sections: %s, entry: %d", fname, nentry);
    }
    return;
}
static void ReadBoundaryData(FILE *fp, Space *space, const int n)
{
    Partition *const part = &(space->part);
    String str = {'\0'}; /* store the current read line */
    const char *fmtI = ParseFormat("%lg");
    ParseCommand(fgets(str, sizeof str, fp));
    if (0 == strncmp(str, "inflow", sizeof str)) {
        part->typeBC[n] = INFLOW;
        ReadConsecutiveData(fp, VARBC - 1, fmtI, part->varBC[n], NULL);
        return;
    }
    if (0 == strncmp(str, "outflow", sizeof str)) {
        part->typeBC[n] = OUTFLOW;
        return;
    }
    if (0 == strncmp(str, "slip wall", sizeof str)) {
        part->typeBC[n] = SLIPWALL;
        Sread(fp, 1, fmtI, &(part->varBC[n][VARBC-1]));
        return;
    }
    if (0 == strncmp(str, "noslip wall", sizeof str)) {
        part->typeBC[n] = NOSLIPWALL;
        Sread(fp, 1, fmtI, &(part->varBC[n][VARBC-1]));
        return;
    }
    if (0 == strncmp(str, "periodic", sizeof str)) {
        part->typeBC[n] = PERIODIC;
        return;
    }
    ShowError("unidentified boundary type: n: %d, type: %s", n, str);
    return;
}
static void ReadConsecutiveData(FILE *fp, const int n, const char *fmt,
        Real *preal, char pstr[][VARSTR])
{
    if (NULL != preal) {
        for (int m = 0; m < n; ++m) {
            Sread(fp, 1, fmt, preal + m);
        }
    } else {
        String str = {'\0'};
        for (int m = 0; m < n; ++m) {
            ParseCommand(fgets(str, sizeof str, fp));
            strncpy(pstr[m], str, sizeof pstr[m]);
        }
    }
    return;
}
static void WriteBoundaryData(FILE *fp, const Space *space, const int n)
{
    const Partition *const part = &(space->part);
    switch (part->typeBC[n]) {
        case INFLOW:
            fprintf(fp, "boundary type: inflow\n");
            fprintf(fp, "density: %.6g\n", part->varBC[n][0]);
            fprintf(fp, "x velocity: %.6g\n", part->varBC[n][1]);
            fprintf(fp, "y velocity: %.6g\n", part->varBC[n][2]);
            fprintf(fp, "z velocity: %.6g\n", part->varBC[n][3]);
            fprintf(fp, "pressure: %.6g\n", part->varBC[n][4]);
            break;
        case OUTFLOW:
            fprintf(fp, "boundary type: outflow\n");
            break;
        case SLIPWALL:
            fprintf(fp, "boundary type: slip wall\n");
            fprintf(fp, "temperature: %.6g\n", part->varBC[n][VARBC-1]);
            break;
        case NOSLIPWALL:
            fprintf(fp, "boundary type: noslip wall\n");
            fprintf(fp, "temperature: %.6g\n", part->varBC[n][VARBC-1]);
            break;
        case PERIODIC:
            fprintf(fp, "boundary type: periodic\n");
            break;
        default:
            ShowError("unidentified boundary type: n: %d, type: %d", n, part->typeBC[n]);
            break;
    }
    return;
}
static void WriteInitializerData(FILE *fp, const Space *space, const int n)
{
    const Partition *const part = &(space->part);
    switch (part->typeIC[n]) {
        case ICGLOBAL:
            break;
        case ICPLANE:
            fprintf(fp, "regional initialization: plane\n");
            fprintf(fp, "plane point x, y, z: %.6g, %.6g, %.6g\n",
                    part->posIC[n][0], part->posIC[n][1], part->posIC[n][2]);
            fprintf(fp, "plane normal nx, ny, nz: %.6g, %.6g, %.6g\n",
                    part->posIC[n][3], part->posIC[n][4], part->posIC[n][5]);
            break;
        case ICSPHERE:
            fprintf(fp, "regional initialization: sphere\n");
            fprintf(fp, "center point x, y, z: %.6g, %.6g, %.6g\n",
                    part->posIC[n][0], part->posIC[n][1], part->posIC[n][2]);
            fprintf(fp, "radius: %.6g\n", part->posIC[n][6]);
            break;
        case ICBOX:
            fprintf(fp, "regional initialization: box\n");
            fprintf(fp, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n",
                    part->posIC[n][0], part->posIC[n][1], part->posIC[n][2]);
            fprintf(fp, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n",
                    part->posIC[n][3], part->posIC[n][4], part->posIC[n][5]);
            break;
        case ICCYLINDER:
            fprintf(fp, "regional initialization: cylinder\n");
            fprintf(fp, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n",
                    part->posIC[n][0], part->posIC[n][1], part->posIC[n][2]);
            fprintf(fp, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n",
                    part->posIC[n][3], part->posIC[n][4], part->posIC[n][5]);
            fprintf(fp, "radius: %.6g\n", part->posIC[n][6]);
            break;
        default:
            break;
    }
    fprintf(fp, "density: %s\n", part->varIC[n][0]);
    fprintf(fp, "x velocity: %s\n", part->varIC[n][1]);
    fprintf(fp, "y velocity: %s\n", part->varIC[n][2]);
    fprintf(fp, "z velocity: %s\n", part->varIC[n][3]);
    fprintf(fp, "pressure: %s\n", part->varIC[n][4]);
    return;
}
static void WriteVerifyData(const Time *time, const Space *space, const Model *model)
{
    const Partition *const part = &(space->part);
    const char *fname = "artracfd.verify";
    FILE *fp = Fopen(fname, "w");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "#                     Case Verification for ArtraCFD                          -\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                          >> Space Domain <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n",
            part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]);
    fprintf(fp, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n",
            part->domain[X][MAX], part->domain[Y][MAX], part->domain[Z][MAX]);
    fprintf(fp, "mx, my, mz: %d, %d, %d\n", part->m[X], part->m[Y], part->m[Z]);
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                          >> Time Domain <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "restart number tag: %d\n", time->restart);
    fprintf(fp, "termination time: %.6g\n", time->end);
    fprintf(fp, "CFL condition number: %.6g\n", time->numCFL);
    fprintf(fp, "maximum computing steps: %d\n", time->stepN);
    fprintf(fp, "space data writing frequency: %d\n", time->dataW[PROSD]);
    fprintf(fp, "data streamer: %d\n", time->dataStreamer);
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                        >> Numerical Method <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "temporal scheme: %d\n", model->tScheme);
    fprintf(fp, "spatial scheme: %d\n", model->sScheme);
    fprintf(fp, "dimensional scheme: %d\n", model->multidim);
    fprintf(fp, "Jacobian average: %d\n", model->jacobMean);
    fprintf(fp, "flux splitting method: %d\n", model->fluxSplit);
    fprintf(fp, "phase interaction: %d\n", model->psi);
    fprintf(fp, "ibm reconstruction layers: %d\n", model->ibmLayer);
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                       >> Material Properties <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "material: %d\n", model->mid);
    fprintf(fp, "viscous level: %.6g\n", model->refMu);
    fprintf(fp, "gravity state: %d\n", model->gState);
    fprintf(fp, "gravity vector: %.6g, %.6g, %.6g\n", model->g[X], model->g[Y], model->g[Z]);
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                        >> Reference Values  <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "length: %.6g\n", model->refL);
    fprintf(fp, "density: %.6g\n", model->refRho);
    fprintf(fp, "velocity: %.6g\n", model->refV);
    fprintf(fp, "temperature: %.6g\n", model->refT);
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                     >> Initialization <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    WriteInitializerData(fp, space, ICGLOBAL);
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                     >> Boundary Condition <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "Domian West\n");
    WriteBoundaryData(fp, space, PWB);
    fprintf(fp, "#\n");
    fprintf(fp, "Domian East\n");
    WriteBoundaryData(fp, space, PEB);
    fprintf(fp, "#\n");
    fprintf(fp, "Domian South\n");
    WriteBoundaryData(fp, space, PSB);
    fprintf(fp, "#\n");
    fprintf(fp, "Domian North\n");
    WriteBoundaryData(fp, space, PNB);
    fprintf(fp, "#\n");
    fprintf(fp, "Domian Front\n");
    WriteBoundaryData(fp, space, PFB);
    fprintf(fp, "#\n");
    fprintf(fp, "Domian Back\n");
    WriteBoundaryData(fp, space, PBB);
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                  >> Regional Initialization <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    for (int n = 1; n < part->nIC; ++n) {
        fprintf(fp, "#\n");
        WriteInitializerData(fp, space, n);
    }
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#                    >> Field Data Probes <<\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "point probe count: %d\n", time->dataN[PROPT]);
    fprintf(fp, "line probe count: %d\n", time->dataN[PROLN]);
    fprintf(fp, "curve probe count: %d\n", time->dataN[PROCV]);
    fprintf(fp, "force probe count: %d\n", time->dataN[PROFC]);
    fprintf(fp, "#\n");
    fprintf(fp, "point probe writing frequency: %d\n", time->dataW[PROPT]);
    fprintf(fp, "line probe writing frequency: %d\n", time->dataW[PROLN]);
    fprintf(fp, "body-conformal probe writing frequency: %d\n", time->dataW[PROCV]);
    fprintf(fp, "surface force writing frequency: %d\n", time->dataW[PROFC]);
    fprintf(fp, "#\n");
    for (int n = 0; n < time->dataN[PROPT]; ++n) {
        fprintf(fp, "point probe x, y, z: %.6g, %.6g, %.6g\n",
                time->pp[n][0], time->pp[n][1], time->pp[n][2]);
    }
    fprintf(fp, "#\n");
    for (int n = 0; n < time->dataN[PROLN]; ++n) {
        fprintf(fp, "line probe x1, y1, z1: %.6g, %.6g, %.6g\n",
                time->lp[n][0], time->lp[n][1], time->lp[n][2]);
        fprintf(fp, "line probe x2, y2, z2: %.6g, %.6g, %.6g\n",
                time->lp[n][3], time->lp[n][4], time->lp[n][5]);
        fprintf(fp, "resolution: %.6g\n", time->lp[n][6]);
    }
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fclose(fp);
    return;
}
static void CheckCaseSettingData(const Time *time, const Space *space, const Model *model)
{
    const Partition *const part = &(space->part);
    const Real zero = 0.0;
    /* space */
    if ((zero >= (part->domain[X][MAX] - part->domain[X][MIN])) ||
            (zero >= (part->domain[Y][MAX] - part->domain[Y][MIN])) ||
            (zero >= (part->domain[Z][MAX] - part->domain[Z][MIN]))) {
        ShowError("domain region should have max > min");
    }
    if ((1 > part->m[X]) || (1 > part->m[Y]) || (1 > part->m[Z])) {
        ShowError("mesh number should be positive");
    }
    if ((1 > part->proc[X]) || (1 > part->proc[Y]) || (1 > part->proc[Z])) {
        ShowError("processor number should be positive");
    }
    /* time */
    if ((0 > time->restart) || (zero >= time->end) || (zero >= time->numCFL)) {
        ShowError("values in time section should not be negative");
    }
    /* numerical method */
    if ((0 > model->tScheme) || (0 > model->sScheme) || (0 > model->multidim) ||
            (0 > model->jacobMean) || (0 > model->fluxSplit) || (0 > model->psi)) {
        ShowError("values in numerical section should not be negative");
    }
    /* material */
    if ((0 > model->mid)) {
        ShowError("material type should not be negative");
    }
    /* reference */
    if ((zero >= model->refL) || (zero >= model->refRho) ||
            (zero >= model->refV) || (zero >= model->refT)) {
        ShowError("reference values should be positive");
    }
    return;
}
/* a good practice: end file with a newline */

