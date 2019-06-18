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

VolumeStruct* ByteVolume(VolumeStruct *in)
{
  VolumeStruct *volume;
  int i;
  ushort max;
  
  volume = calloc(1, sizeof(VolumeStruct));

  volume->VolWidth  = in->VolWidth;
  volume->VolHeight = in->VolHeight;
  volume->VolDepth  =in->VolDepth;
  volume->VolSize   = in->VolSize;
  volume->DataType=1;

  volume->data.uchar_data = (uchar*)calloc(volume->VolSize, sizeof(uchar));
  max=0;
  for (i=0;i<in->VolSize;i++)
    if (max<in->data.ushort_data[i])
      max=in->data.ushort_data[i];
  
  for (i=0;i<in->VolSize;i++)
    volume->data.uchar_data[i]=(uchar)(255.0*(float)(in->data.ushort_data[i])/((float)max));

  return volume; 

}

int main(int argc, char **argv)
{
  VolumeStruct VolStruct, *out;
  avs_header avs_head;
  int endian_type = 0;
    
  if (argc < 2) {
    printf("Usage: %s <infile.fld> <outfile.fld> \n", argv[0]);
    exit(1);
  }
  
  
  if(!readAVS(argv[1], &avs_head, &VolStruct)) 
     return(0);
 
  if (avs_head.datatype>1)
    {
      out =  ByteVolume(&VolStruct);
      avs_head.datatype=1;
      writeAVS(argv[2], &avs_head, out);
    }
  return(1);
}
 
