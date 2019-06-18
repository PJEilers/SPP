/*************************************************************************************/
/* ReadWriteData.h - Header file for ReadWriteData.c                                 */
/* Usage     :                                                                       */
/* Author    :   Georgios K. Ouzounis                                                */
/* Dep.      :   IWI - University of Groningen                                       */
/* Vs&Date   :   vs 1.0 - 10th Feb. 2005                                             */  
/*************************************************************************************/

#ifndef _READWRITEDATA_H_
#define _READWRITEDATA_H_

#include "types.h"

#ifdef _cplusplus
extern "C" {
#endif

void readByteData(FILE *f, uchar *data, int xdim, int ydim, int zdim);
void readShortData(FILE *f, ushort *data, int xdim, int ydim, int zdim);
void writeByteData(FILE *f, uchar *data, int xdim, int ydim, int zdim);
void writeShortData(FILE *f, ushort *data, int xdim, int ydim, int zdim);
#ifdef _cplusplus
}
#endif

#endif
