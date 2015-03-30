/****************************************************************************
 * Export and Load Computed Data                                            *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "datastream.h"
#include <stdio.h> /* standard library for input and output */
#include "paraview.h"
#include "parasight.h"
#include "ensight.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteComputedData(const Real *U, const Space *space, const Time *time, 
        const Partition *part, const Flow *flow)
{
    switch (time->dataStreamer) {
        case 0: /* ParaView */
            WriteComputedDataParaview(U, space, time, part, flow);
            break;
        case 1: /* Generic Ensight */
            WriteComputedDataEnsight(U, space, time, part, flow);
            break;
        case 2: /* ParaView Ensight */
            WriteComputedDataParasight(U, space, time, part, flow);
            break;
        default: /* ParaView */
            WriteComputedDataParaview(U, space, time, part, flow);
            break;
    }
    return 0;
}
int LoadComputedData(Real *U, const Space *space, Time *time,
        const Partition *part, const Flow *flow)
{
    switch (time->dataStreamer) {
        case 0: /* ParaView */
            LoadComputedDataParaview(U, space, time, part, flow);
            break;
        case 1: /* Ensight */
            LoadComputedDataEnsight(U, space, time, part, flow);
            break;
        case 2: /* Paraview Ensight */
            LoadComputedDataParasight(U, space, time, part, flow);
            break;
        default: /* ParaView */
            LoadComputedDataParaview(U, space, time, part, flow);
            break;
    }
    return 0;
}
/* a good practice: end file with a newline */

