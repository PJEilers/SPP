/************************************************************************************/
/* sff_io.h  - Header file for sff_io.c                                             */
/* Usage   :                                                                        */
/* Author  :   Georgios K. Ouzounis                                                 */
/* Dep.    :   IWI - University of Groningen                                        */
/* Vs&Date :   vs 1.0 - 2nd Nov. 2004                                               */  
/************************************************************************************/

#ifndef _SFF_IO_H_
#define _SFF_IO_H_

#include "types.h"
#include "Volume.h"


typedef struct sff_header
{
   uchar   datatype;            /* short or byte  */
   uchar   ndim;                /* number of dims */
   int     dim1, dim2, dim3;    /* volume dims    */               
} sff_header;

void sff_read_header(FILE *f, sff_header *header);
int readSFF(char *in_file, sff_header *header, VolumeStruct *VolStruct, int endian_type);
void sff_write_header(FILE *f, sff_header *header); 
int writeSFF(char *out_file, sff_header *header, VolumeStruct *VolStruct);


#endif
