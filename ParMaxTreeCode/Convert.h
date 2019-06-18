/*************************************************************************************/
/* convert.h :   Header file for Convert.c                                           */
/* Usage     :                                                                       */
/* Author    :   Georgios K. Ouzounis                                                */
/* Dep.      :   IWI - University of Groningen                                       */
/* Vs&Date   :   vs 1.0 - 10th Feb. 2005                                             */  
/*************************************************************************************/

#ifndef _CONVERT_H_
#define _CONVERT_H_

#include "Volume.h"

#ifdef _cplusplus
extern "C" {
#endif

void ConvertByte2Short(VolumeStruct *volume);
void ConvertShort2Byte(VolumeStruct *volume);

#ifdef _cplusplus
}
#endif

#endif
