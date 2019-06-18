/*******************************************************************/
/* RAW2AVS.c - Performs volume convertion from .sff to AVS .fld    */
/* Usage   :   sff2avs infile.raw xdim ydim zdim type outfile.fld  */
/* Specs.  :   Reads unsigned char and unsigned short datatypes    */
/*             Exports the same datatype.                          */
/*             Supports uniform fields with 1 component/element    */
/*             i.e. veclen =1 and produces a single file output    */
/* Author  :   Michael H.F. Wilkinson                              */
/* Dep.    :   IWI - University of Groningen                       */
/* Vs&Date :   vs 1.0 - 28th Oct. 2004                             */  
/*******************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include "Volume.h"
#include "ReadWriteData.h"
#include "avs_io.h"

int main(int argc, char **argv)
{
  VolumeStruct *vol_struct=malloc(sizeof(VolumeStruct));
  int xdim, ydim, zdim,dtype,i;
  avs_header avs_head;
  FILE *f;
  uchar temp, *pchar;
  
  int endian_type=0;
  int byteoffset=0;    
  if (argc < 7) {
    printf("Usage: %s <infile.raw> <xsize> <ysize> <zsize> <dtype> <outfile.fld>[byte offset]\n", argv[0]);
    exit(1);
  }

  f = fopen(argv[1], "r");
  
  if(!f)
  {
      printf("raw2avs: Cannot open the source file\n");
      return(0);
  }
  
  xdim = atoi(argv[2]);
  ydim = atoi(argv[3]);
  zdim = atoi(argv[4]);
  dtype = atoi(argv[5]);
  if (argc > 7){
    byteoffset=atoi(argv[7]);
    fseek(f,byteoffset,SEEK_SET);
  }

  InitVolume(vol_struct, xdim, ydim, zdim, dtype);

  /*if(!readSFF(argv[1], &sff_head, &VolStruct, endian_type)) 
    return(0);
  */
  switch (dtype) {
  case 1: /* bytes */
    readByteData(f, vol_struct->data.uchar_data, xdim, ydim, zdim);
    printf("raw2avs: DataType: Bytes\n");
    break;    
  case 2: /* shorts */    
  case 3: /* shorts */    
    readShortData(f, vol_struct->data.ushort_data, xdim, ydim,zdim);
    /* if Big Endian, swap byte order */
    if(endian_type!=0)
    {
       for(i=0;i<vol_struct->VolSize;i++)
       {
           
	   pchar = (uchar*)&vol_struct->data.ushort_data[i];
	   temp = pchar[0];
	   pchar[0] = pchar[1];
	   pchar[1] = temp;
       }
    }
    printf("raw2avs: DataType: Shorts\n");
    break;     
  default:
    printf("raw2avs: Wrong datatype\n");
    DestroyVolume(vol_struct);
    fclose(f);
    return 0;
  }

 
  /* Prepare avs_header fields */ 
  avs_head.ndim = 3;   
  avs_head.dim1 = xdim;
  avs_head.dim2 = ydim;
  avs_head.dim3 = zdim;
  avs_head.min_x = 0;
  avs_head.min_y = 0;
  avs_head.min_z = 0;
  avs_head.max_x = avs_head.dim1-1;
  avs_head.max_y = avs_head.dim2-1;
  avs_head.max_z = avs_head.dim3-1;
  avs_head.datatype = dtype;
  avs_head.filetype = 0;
  avs_head.skip = 0;
  avs_head.nspace = avs_head.ndim;
  avs_head.veclen = 1;
  avs_head.dataname[0] = '\0';
  
  writeAVS(argv[6], &avs_head, vol_struct);
  
  return(1);
}
 
