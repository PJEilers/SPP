/************************************************************************************/
/* sff_io.c  - Reads and writes SFF volume files (.sff)                             */
/* Usage   :   Implements Read/Write functions for unsigned char and shorts         */
/* Read    :   endian_type: if data are shorts use 0 for little endian, 1 otherwise */
/* Credits :   Always returns little endian files                                   */
/* Author  :   Georgios K. Ouzounis                                                 */
/* Dep.    :   IWI - University of Groningen                                        */
/* Vs&Date :   vs 1.0 - 2nd Nov. 2004                                               */  
/************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "ReadWriteData.h"
#include "sff_io.h"
#include "Volume.h"


/*****************************************************************/
/*                      File IO                                  */
/*****************************************************************/

unsigned int to_int(uchar bytes[4])
{
  return bytes[3]+256*(bytes[2]+256*(bytes[1]+256*bytes[0]));
}
 
/*****************************************************************/
/* AVS Field WRITE Routines                                      */
/*****************************************************************/

void sff_write_header(FILE *f, sff_header *header)  
{ 
  uchar delim;
  fwrite(header, sizeof(*header), 1, f);
  delim = '\f';
  fwrite(&delim, sizeof(uchar), 1, f);
}

int writeSFF(char *out_file, sff_header *header, VolumeStruct *VolStruct)
{
  FILE *f;

  f = fopen(out_file, "w");
  if(!f)
  {
     printf("Cannot open file: %s for writing the AVS field\n", out_file);
     return(0);
  }  
  switch(header->datatype){
     case 1:
       header->datatype = 1;
       sff_write_header(f, header);
       writeByteData(f, VolStruct->data.uchar_data, header->dim1, header->dim2, header->dim3);
       break;
     case 2:
     case 3:
       header->datatype = 3;
       sff_write_header(f, header);
       writeShortData(f, VolStruct->data.ushort_data, header->dim1, header->dim2, header->dim3);       
       break;
     default:
       printf("writeAVS: Datatype is not supported\n");
       fclose(f);
       return(0);
   } 
  fclose(f);
  return(1);
}

/*****************************************************************/
/* SFF READ Routines                                             */
/*****************************************************************/

void read_sff_header(FILE *f, sff_header *header)
{
  unsigned char ndim, ftype, bytes[4], delim;
    
  fread(&ftype, sizeof(unsigned char), 1, f);
  header->datatype = ftype;
  fread(&ndim, sizeof(unsigned char), 1, f);
  if (ndim != 3) 
  {
    printf("The .sff file is not a volume !\n");
    return;
  }
  header->ndim = ndim;

  fread(bytes, sizeof(unsigned char), 4, f);
  header->dim1 = to_int(bytes);
  fread(bytes, sizeof(unsigned char), 4, f);
  header->dim2 = to_int(bytes);
  fread(bytes, sizeof(unsigned char), 4, f);
  header->dim3 = to_int(bytes); 

  /* skip comments */
  delim = fgetc(f);
  while (delim != '\f')
    delim = fgetc(f); 
}

int readSFF(char *filename, sff_header *header, VolumeStruct *vol_struct, int endian_type)
{
  FILE *f;
  uchar temp, *pchar;
  int i;
   
  f = fopen(filename, "r");
  if(!f)
  {
      printf("readSFF: Cannot open the source file\n");
      return(0);
  }
  read_sff_header(f, header);

  InitVolume(vol_struct, header->dim1, header->dim2, header->dim3, header->datatype);
   
  switch (header->datatype) {
  case 1: /* bytes */
    readByteData(f, vol_struct->data.uchar_data, header->dim1, header->dim2, header->dim3);
    printf("readSFF: DataType: Bytes\n");
    break;    
  case 2: /* shorts */    
  case 3: /* shorts */    
    readShortData(f, vol_struct->data.ushort_data, header->dim1, header->dim2, header->dim3);
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
    printf("readSFF: DataType: Shorts\n");
    break;     
  default:
    printf("readSFF: Wrong datatype\n");
    DestroyVolume(vol_struct);
    fclose(f);
    return 0;
  }
  fclose(f); 
  return(1);
}


