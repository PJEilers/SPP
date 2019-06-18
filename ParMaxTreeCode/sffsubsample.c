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

int main(int argc, char **argv)
{
  VolumeStruct VolStruct;
  VolumeStruct *SubSample;

  sff_header sff_head;
  int endian_type = 1;
  int x,y,z,posS,posV;
 
  
  if (argc < 2) {
    printf("Usage: %s <infile.sff> <outfile.sff> \n", argv[0]);
    exit(1);
  }
  
  printf("%s\n",argv[1]);
  
  if(!readSFF(argv[1], &sff_head, &VolStruct, endian_type)) 
     return(0);

  printf("creating volume\n");
  SubSample =  CreateVolume(VolStruct.VolWidth / 2, 
			    VolStruct.VolHeight / 2, 
			    VolStruct.VolDepth / 2, 
			    VolStruct.DataType);

  switch(sff_head.datatype){
  case 1:
    posS=0;
    posV=0;
    for (z = 0; z < VolStruct.VolDepth; z+=2, 
	   posV+= (VolStruct.VolWidth*VolStruct.VolHeight))
      for (y = 0; y < VolStruct.VolHeight; y+=2, posV+= (VolStruct.VolWidth))
	for (x = 0; x < VolStruct.VolWidth; x+=2, posV+=2, posS++){
          /*printf("%d,%d,%d",x,y,z);*/        
	  SubSample->data.uchar_data[posS]=
	    VolStruct.data.uchar_data[posV]; 
	}
    break;
  case 2:
  case 3:
    posS=0;
    posV=0;
    for (z = 0; z < VolStruct.VolDepth; z+=2, 
	   posV+= (VolStruct.VolWidth*VolStruct.VolHeight))
      for (y = 0; y < VolStruct.VolHeight; y+=2, posV+= (VolStruct.VolWidth))
	for (x = 0; x < VolStruct.VolWidth; x+=2, posV+=2, posS++){
          /*printf("%d,%d,%d",x,y,z);  */      
	  SubSample->data.ushort_data[posS]=
	    VolStruct.data.ushort_data[posV]; 
	}
    break;
  default: 
    printf("unknown data type\n");
    exit(1);
  }

  sff_head.dim1 /=2;  
  sff_head.dim2 /=2;  
  sff_head.dim3 /=2;  
 
  writeSFF(argv[2], &sff_head, SubSample);
  
  
  
  return(1);
}
 
