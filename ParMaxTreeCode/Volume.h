/*************************************************************************************/
/* Volume.h :   Header file for Volume.c                                             */
/* Usage    :                                                                        */
/* Author   :   Georgios K. Ouzounis                                                 */
/* Dep.     :   IWI - University of Groningen                                        */
/* Vs&Date  :   vs 1.0 - 14th Jan. 2005                                              */  
/*************************************************************************************/

#ifndef _VOLUME_H_
#define _VOLUME_H_

#include "types.h"

#ifdef _cplusplus
extern "C" {
#endif

typedef union {
    ushort *ushort_data;
    uchar *uchar_data;
} OutData;

typedef struct VolumeStruct {
  long   VolWidth;          /* dim X                */
  long   VolHeight;         /* dim Y                */
  long   VolDepth;          /* dim Z                */
  long   VolSize;           /* Volume Size          */
  int    DataType;          /* either short or byte */

  OutData data;
  
} VolumeStruct;


void InitVolume(VolumeStruct *vol, int width, int height, int depth, int datatype);
VolumeStruct* CreateVolume(int width, int height, int depth, int datatype);
void DeleteVolume(VolumeStruct *volume);
void DestroyVolume(VolumeStruct *volume);

#ifdef _cplusplus
}
#endif

#endif
