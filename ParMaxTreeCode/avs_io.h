/************************************************************************************/
/* avs_io.h  - Header file for avs_io.c                                             */
/* Usage   :                                                                        */
/* Author  :   Georgios K. Ouzounis                                                 */
/* Dep.    :   IWI - University of Groningen                                        */
/* Vs&Date :   vs 1.0 - 2nd Nov. 2004                                               */  
/************************************************************************************/

#ifndef _AVS_IO_H_
#define _AVS_IO_H_

#include "types.h"
#include "Volume.h"


typedef struct avs_header
{
   int   ndim;
   int   dim1, dim2, dim3;    /* volume dims */
   float min_x,min_y,min_z;   /* minmum extend  */
   float max_x,max_y,max_z;   /* maximum extend */
   int   datatype;            /* either short or byte - from SFF file */
   int   filetype;            /* 0: binary, 1: ascii */  
   int   skip;                /* not used here */
   int   nspace;
   int   veclen;
   char  dataname[512];       /* file containing data, NULL otherwise */
} avs_header;


int readAVS(char *in_file, avs_header *header, VolumeStruct *VolStruct);
int writeAVS(char *out_file, avs_header *header, VolumeStruct *VolStruct);

#endif
