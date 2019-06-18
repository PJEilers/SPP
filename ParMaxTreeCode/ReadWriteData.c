/*************************************************************************************/
/* ReadWriteData.c - Reads/Writes from and to 1D uchar or ushort arrays              */
/* Usage     :                                                                       */
/* Author    :   Georgios K. Ouzounis                                                */
/* Dep.      :   IWI - University of Groningen                                       */
/* Vs&Date   :   vs 1.0 - 10th Feb. 2005                                             */  
/*************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "ReadWriteData.h"

void readByteData(FILE *f, uchar *data, int xdim, int ydim, int zdim)
{
  int z, no_read=0, cur_read;

  for (z=0; z<zdim; z++) {
    cur_read = fread(&data[no_read], sizeof(uchar), xdim*ydim, f);
    assert(cur_read==xdim*ydim);
    no_read+=cur_read;
  }
}
void readShortData(FILE *f, ushort *data, int xdim, int ydim, int zdim)
{
  int z, no_read=0, cur_read;
  
  for (z=0; z<zdim; z++) {
    cur_read = fread(&data[no_read], sizeof(ushort), xdim*ydim, f);    
    assert(cur_read==xdim*ydim); 
    no_read+=cur_read;
  }
}
void writeByteData(FILE *f, uchar *data, int xdim, int ydim, int zdim)
{
  int z, no_read=0, cur_read;

  for (z=0; z<zdim; z++) {
    cur_read = fwrite(&data[no_read], sizeof(uchar), xdim*ydim, f);
    assert(cur_read==xdim*ydim);
    no_read+=cur_read;
  }
}
void writeShortData(FILE *f, ushort *data, int xdim, int ydim, int zdim)
{
  int z, no_read=0, cur_read;

  for (z=0; z<zdim; z++) {
    cur_read = fwrite(&data[no_read], sizeof(ushort), xdim*ydim, f);    
    assert(cur_read==xdim*ydim);
    no_read+=cur_read;
  }

}

