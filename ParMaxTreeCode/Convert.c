/*************************************************************************************/
/* convert.c :   Converts byte arrays to short and the opposite                      */
/* Usage     :   Used with volume structures - see: Volume.c                         */
/* Author    :   Georgios K. Ouzounis                                                */
/* Dep.      :   IWI - University of Groningen                                       */
/* Vs&Date   :   vs 1.0 - 10th Feb. 2005                                             */  
/*************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Volume.h"
#include "types.h"

void ConvertByte2Short(VolumeStruct *volume)
{
    int i;
    ushort *ushort_data;
    
    ushort_data = (ushort*)calloc(volume->VolSize, sizeof(ushort));
    for(i=0;i<volume->VolSize;i++)
       ushort_data[i] = (ushort)volume->data.uchar_data[i];
    volume->DataType = 2;
    free(volume->data.uchar_data);
    volume->data.ushort_data = ushort_data;           
}

void ConvertShort2Byte(VolumeStruct *volume)
{
    int i;
    uchar *uchar_data;
    
    uchar_data = (uchar*)calloc(volume->VolSize, sizeof(uchar));
    for(i=0;i<volume->VolSize;i++)
       uchar_data[i] = (uchar)volume->data.ushort_data[i];
    volume->DataType = 1;
    free(volume->data.ushort_data);
    volume->data.uchar_data = uchar_data;           
}

