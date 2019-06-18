/*****************************************************************/
/* SFF2AVS.c - Performs volume convertion from .sff to AVS .fld  */
/* Usage   :   sff2avs infile.sff outfile.fld                    */
/* Specs.  :   Reads unsigned char and unsigned short datatypes  */
/*             Exports the same datatype.                        */
/*             Supports uniform fields with 1 component/element  */
/*             i.e. veclen =1 and produces a single file output  */
/* Author  :   Georgios K. Ouzounis                              */
/* Dep.    :   IWI - University of Groningen                     */
/* Vs&Date :   vs 1.0 - 28th Oct. 2004                           */  
/*****************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include "Volume.h"
#include "sff_io.h"
#include "avs_io.h"

int main(int argc, char **argv)
{
  VolumeStruct VolStruct;
  avs_header avs_head;
  sff_header sff_head;
  int endian_type;
    
  if (argc < 2) {
    printf("Usage: %s <infile.sff> <outfile.fld> <0: little endian OR 1: big endian>\n", argv[0]);
    exit(1);
  }
  
  endian_type = atoi((const char*)argv[3]);
  
  if(!readSFF(argv[1], &sff_head, &VolStruct, endian_type)) 
     return(0);
 
  /* Prepare avs_header fields */ 
  avs_head.ndim = sff_head.ndim;   
  avs_head.dim1 = sff_head.dim1;
  avs_head.dim2 = sff_head.dim2;
  avs_head.dim3 = sff_head.dim3;
  avs_head.min_x = 0;
  avs_head.min_y = 0;
  avs_head.min_z = 0;
  avs_head.max_x = avs_head.dim1-1;
  avs_head.max_y = avs_head.dim2-1;
  avs_head.max_z = avs_head.dim3-1;
  avs_head.datatype = sff_head.datatype;
  avs_head.filetype = 0;
  avs_head.skip = 0;
  avs_head.nspace = avs_head.ndim;
  avs_head.veclen = 1;
  avs_head.dataname[0] = '\0';
  
  writeAVS(argv[2], &avs_head, &VolStruct);
  
  return(1);
}
 
