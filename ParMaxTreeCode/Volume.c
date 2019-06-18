/*************************************************************************************/
/* Volume.c :   Volume manipulation routines                                         */
/* Usage    :                                                                        */
/* Author   :   Georgios K. Ouzounis                                                 */
/* Dep.     :   IWI - University of Groningen                                        */
/* Vs&Date  :   vs 1.0 - 14th Jan. 2005                                              */  
/*************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Volume.h"

void InitVolume(VolumeStruct *volume, int width, int height, int depth, int datatype)
{
  volume->VolWidth  = width;
  volume->VolHeight = height;
  volume->VolDepth  = depth;
  volume->VolSize   = width*height*depth;
  volume->DataType  = datatype;

  if(volume->DataType==1) /* uchar */
    volume->data.uchar_data = (uchar*)calloc(volume->VolSize, sizeof(uchar));
  else /* ushort */
    volume->data.ushort_data = (ushort*)calloc(volume->VolSize, sizeof(ushort));
}

VolumeStruct* CreateVolume(int width, int height, int depth, int datatype)
{
  VolumeStruct *volume;
  
  volume = calloc(1, sizeof(VolumeStruct));

  volume->VolWidth  = width;
  volume->VolHeight = height;
  volume->VolDepth  = depth;
  volume->VolSize   = width*height*depth;
  volume->DataType  = datatype;

  if(volume->DataType==1) /* uchar */
    volume->data.uchar_data = (uchar*)calloc(volume->VolSize, sizeof(uchar));
  else /* ushort */
    volume->data.ushort_data = (ushort*)calloc(volume->VolSize, sizeof(ushort));
  
  return volume; 

}
/*************************************************************************************/
/* If volume is not a pointer (the structure is kept on the stack) use this function */
/*************************************************************************************/
void DestroyVolume(VolumeStruct *volume)
{
  if(volume) {
    if(volume->DataType==1) /* uchar */
      free(volume->data.uchar_data);
    else /* ushort */
      free(volume->data.ushort_data);
    
    volume=0;
  }
}

void DeleteVolume(VolumeStruct *volume)
{
  if(volume) {
    if(volume->DataType==1) /* uchar */
      free(volume->data.uchar_data);
    else /* ushort */
      free(volume->data.ushort_data);
    free(volume);
    volume=0;
  }
}
