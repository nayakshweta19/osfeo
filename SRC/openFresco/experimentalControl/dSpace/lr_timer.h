/* lr_timer.h ***********************************************************
*                                                                       *
* Low-Resolution Timer Utilities                                        *
*                                                                       *
* uses the PC/AT's real-time clock to calibrate software delay          *
* 1 msec resolution                                                     *
*                                                                       *
* 1993 - 1996 by dSPACE GmbH, Paderborn                                        *
************************************************************************/

/* $RCSfile$ $Revision: 314 $ $Date: 2011-05-23 05:17:07 +0800 (星期一, 23 五月 2011) $ */

#ifndef _LR_TIMER
#define _LR_TIMER

#ifdef __cplusplus
extern "C" {
#endif

#include <dstypes.h>

void LR_TIMER_delay (UInt32 duration);             /* duration in msec */

#ifdef __cplusplus
}
#endif

#endif
