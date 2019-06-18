#include <stddef.h>
#include <stdlib.h>
#include <strings.h>
#include "macros.h"
#include "types.h"
#include "Volume.h"
#include "se.h"

void DilateHorLine(ushort *f, ulong width, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r)
/* k is length of SE in number of pixels */
/* width is length of g, h, h2, and r */
{
   ulong x;

   for (x=0; x<width; x++)
   {
      if (x%k)  g[x] = MAX(g[x-1], f[x]);
      else  g[x] = f[x];
   }
   for (x=width-1; x>0; x--)
   {
      if (((x%k)==(k-1)) || (x==(width-1)))  h[x] = f[x];
      else  h[x] = MAX(h[x+1], f[x]);
   }
   if ((k==1) || (width==1))  h[0] = f[0];
   else  h[0] = MAX(h[1], f[0]);
   h2[width-1] = f[width-1];
   for (x=width-2; x>=(width-k); x--)  h2[x] = MAX(h2[x+1], f[x]);
   h2[0] = MAX(h2[1], f[0]);

   if (width <= (k/2))
   {
      for (x=0; x<width; x++)  r[x] = g[width-1];
   }
   else if (width <= k)
   {
      for (x=0; x<(width-k/2); x++)  r[x] = g[x+k/2];
      for (; x<=(k/2); x++)  r[x] = g[width-1];
      for (; x<width; x++)  r[x] = h[x-k/2];
   }
   else  /* width > k */
   {
      for (x=0; x<(width-k/2); x++)
      {
         if (x < (k/2))  r[x] = g[x+k/2];
         else  r[x] = MAX(g[x+k/2],h[x-k/2]);
      }
      for (x=width-k/2; x<width; x++)  r[x] = h2[x-k/2];
   }
} /* DilateHorLine */

void ErodeHorLine(ushort *f, ulong width, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r)
/* k is length of SE in number of pixels */
/* width is length of g, h, h2, and r */
{
   ulong x;

   for (x=0; x<width; x++)
   {
      if (x%k)  g[x] = MIN(g[x-1], f[x]);
      else  g[x] = f[x];
   }
   for (x=width-1; x>0; x--)
   {
      if (((x%k)==(k-1)) || (x==(width-1)))  h[x] = f[x];
      else  h[x] = MIN(h[x+1], f[x]);
   }
   if ((k==1) || (width==1))  h[0] = f[0];
   else  h[0] = MIN(h[1], f[0]);
   h2[width-1] = f[width-1];
   for (x=width-2; x>=(width-k); x--)  h2[x] = MIN(h2[x+1], f[x]);
   h2[0] = MIN(h2[1], f[0]);

   if (width <= (k/2))
   {
      for (x=0; x<width; x++)  r[x] = g[width-1];
   }
   else if (width <= k)
   {
      for (x=0; x<(width-k/2); x++)  r[x] = g[x+k/2];
      for (; x<=(k/2); x++)  r[x] = g[width-1];
      for (; x<width; x++)  r[x] = h[x-k/2];
   }
   else  /* width > k */
   {
      for (x=0; x<(width-k/2); x++)
      {
         if (x < (k/2))  r[x] = g[x+k/2];
         else  r[x] = MIN(g[x+k/2],h[x-k/2]);
      }
      for (x=width-k/2; x<width; x++)  r[x] = h2[x-k/2];
   }
} /* ErodeHorLine */



void DilateVerLine(ushort *f, ulong width, ulong height, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r)
/* k is length of SE in number of pixels */
/* height is length of g, h, h2, and r */
{
   ulong y;

   for (y=0; y<height; y++)
   {
      if(y%k)  
         g[y] = MAX(g[y-1], f[y*width]);
      else  
         g[y] = f[y*width];
   }
   for (y=height-1; y>0; y--)
   {
      if(((y%k)==(k-1)) || (y==(height-1)))  
         h[y] = f[y*width];
      else  
         h[y] = MAX(h[y+1], f[y*width]);
   }
   if ((k==1) || (height==1))  
       h[0] = f[0];
   else  
       h[0] = MAX(h[1], f[0]);
   
   h2[height-1] = f[(height-1)*width];
   
   for(y=height-2; y>=(height-k); y--)  
      h2[y] = MAX(h2[y+1], f[y*width]);
   
   h2[0] = MAX(h[1], f[0]);

   if (height <= (k/2))
   {
      for (y=0; y<height; y++)  r[y*width] = g[height-1];
   }
   else if (height <= k)
   {
      for(y=0; y<(height-k/2); y++)  
          r[y*width] = g[y+k/2];
      for(; y<=(k/2); y++)  
          r[y*width] = g[height-1];
      for(; y<height; y++)  
          r[y*width] = h[y-k/2];
   }
   else  /* height > k */
   {
      for (y=0; y<(height-k/2); y++)
      {
         if(y < (k/2))  
	    r[y*width] = g[y+k/2];
         else
	    r[y*width] = MAX(g[y+k/2],h[y-k/2]);
      }
      for(y=height-k/2; y<height; y++)  
          r[y*width] = h2[y-k/2];
   }
} /* DilateVerLine */



void ErodeVerLine(ushort *f, ulong width, ulong height, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r)
/* k is length of SE in number of pixels */
/* height is length of g, h, and r */
{
   ulong y;

   for (y=0; y<height; y++)
   {
      if (y%k)  g[y] = MIN(g[y-1], f[y*width]);
      else  g[y] = f[y*width];
   }
   for (y=height-1; y>0; y--)
   {
      if (((y%k)==(k-1)) || (y==(height-1)))  h[y] = f[y*width];
      else  h[y] = MIN(h[y+1], f[y*width]);
   }
   if ((k==1) || (height==1))  h[0] = f[0];
   else  h[0] = MIN(h[1], f[0]);
   h2[height-1] = f[(height-1)*width];
   for (y=height-2; y>=(height-k); y--)  h2[y] = MIN(h2[y+1], f[y*width]);
   h2[0] = MIN(h2[1], f[0]);

   if (height <= (k/2))
   {
      for (y=0; y<height; y++)  r[y*width] = g[height-1];
   }
   else if (height <= k)
   {
      for (y=0; y<(height-k/2); y++)  r[y*width] = g[y+k/2];
      for (; y<=(k/2); y++)  r[y*width] = g[height-1];
      for (; y<height; y++)  r[y*width] = h[y-k/2];
   }
   else  /* height > k */
   {
      for (y=0; y<(height-k/2); y++)
      {
         if (y < (k/2))  r[y*width] = g[y+k/2];
         else  r[y*width] = MIN(g[y+k/2],h[y-k/2]);
      }
      for (y=height-k/2; y<height; y++)  r[y*width] = h2[y-k/2];
   }
} /* ErodeVerLine */

void DilateDepthLine(ushort *f, ulong width, ulong height, ulong depth, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r)
/* k is length of SE in number of pixels */
/* depth is length of g, h, h2, and r */
{
   ulong z;

   for (z=0; z<depth; z++)
   {
      if (z%k)  g[z] = MAX(g[z-1], f[width*(height*z)]);
      else  g[z] = f[width*(height*z)];
   }
   for (z=depth-1; z>0; z--)
   {
      if (((z%k)==(k-1)) || (z==(depth-1)))  h[z] = f[width*(height*z)];
      else  h[z] = MAX(h[z+1], f[width*(height*z)]);
   }
   if ((k==1) || (depth==1))  
      h[0] = f[0];
   else  
      h[0] = MAX(h[1], f[0]);
   
   h2[depth-1] = f[(depth-1)*height*width];
   
   for (z=depth-2; z>=(depth-k); z--)  
      h2[z] = MAX(h2[z+1], f[width*(height*z)]);
   h2[0] = MAX(h[1], f[0]);

   if (depth <= (k/2))
   {
      for (z=0; z<depth; z++)  r[width*(height*z)] = g[depth-1];
   }
   else if (depth <= k)
   {
      for (z=0; z<(depth-k/2); z++)  r[width*(height*z)] = g[z+k/2];
      for (; z<=(k/2); z++)  r[width*(height*z)] = g[depth-1];
      for (; z<depth; z++)  r[width*(height*z)] = h[z-k/2];
   }
   else  /* height > k */
   {
      for (z=0; z<(depth-k/2); z++)
      {
         if (z < (k/2))  r[width*(height*z)] = g[z+k/2];
         else  r[width*(height*z)] = MAX(g[z+k/2],h[z-k/2]);
      }
      for (z=depth-k/2; z<depth; z++)  r[width*(height*z)] = h2[z-k/2];
   }
} /* DilateDepthLine */



void ErodeDepthLine(ushort *f, ulong width, ulong height, ulong depth, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r)
/* k is length of SE in number of pixels */
/* depth is length of g, h, and r */
{
   ulong z;

   for (z=0; z<depth; z++)
   {
      if (z%k)  g[z] = MIN(g[z-1], f[width*(height*z)]);
      else  g[z] = f[width*(height*z)];
   }
   for (z=depth-1; z>0; z--)
   {
      if (((z%k)==(k-1)) || (z==(depth-1)))  h[z] = f[width*(height*z)];
      else  h[z] = MIN(h[z+1], f[width*(height*z)]);
   }
   if ((k==1) || (depth==1))  
      h[0] = f[0];
   else  
      h[0] = MIN(h[1], f[0]);
   
   h2[depth-1] = f[(depth-1)*height*width];
   
   for (z=depth-2; z>=(depth-k); z--)  
      h2[z] = MIN(h2[z+1], f[width*(height*z)]);
   h2[0] = MIN(h[1], f[0]);

   if (depth <= (k/2))
   {
      for (z=0; z<depth; z++)  r[width*(height*z)] = g[depth-1];
   }
   else if (depth <= k)
   {
      for (z=0; z<(depth-k/2); z++)  r[width*(height*z)] = g[z+k/2];
      for (; z<=(k/2); z++)  r[width*(height*z)] = g[depth-1];
      for (; z<depth; z++)  r[width*(height*z)] = h[z-k/2];
   }
   else  /* height > k */
   {
      for (z=0; z<(depth-k/2); z++)
      {
         if (z < (k/2))  r[width*(height*z)] = g[z+k/2];
         else  r[width*(height*z)] = MIN(g[z+k/2],h[z-k/2]);
      }
      for (z=depth-k/2; z<depth; z++)  r[width*(height*z)] = h2[z-k/2];
   }
} /* ErodeDepthLine */

void VolGrayDilateHor(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out)
{
    ushort *f = vol_in->data.ushort_data, *r = vol_out->data.ushort_data;
    ulong width = vol_in->VolWidth, depth = vol_in->VolDepth, y, z;
    for(z=0; z<depth; z++)
    {
       for(y=0; y<vol_in->VolHeight; y++)
       {
          DilateHorLine(f, width, k, g, h, h2, r);
          f += width;
          r += width;
       }
    }
} /* ImageGrayDilateHor */

void VolGrayErodeHor(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out)
{
    ushort *f = vol_in->data.ushort_data, *r = vol_out->data.ushort_data;
    ulong width = vol_in->VolWidth, depth = vol_in->VolDepth, y, z;
    for(z=0; z<depth; z++)
    {
       for(y=0; y<vol_in->VolHeight; y++)
       {
          ErodeHorLine(f, width, k, g, h, h2, r);
          f += width;
          r += width;
       }
    }
} /* ImageGrayErodeHor */



void VolGrayDilateVer(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out)
{
    ushort *f = vol_in->data.ushort_data, *r = vol_out->data.ushort_data;
    ulong width = vol_in->VolWidth, height = vol_in->VolHeight, depth = vol_in->VolDepth, x, z;

    for(z=0; z<depth; z++)
    {
       for(x=0; x<width; x++)
       {
          DilateVerLine(f, width, height, k, g, h, h2, r);
          f++;
          r++;
       }
       f += width*(height-1);
       r += width*(height-1);
    }
} /* ImageGrayDilateVer */



void VolGrayErodeVer(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out)
{
    ushort *f = vol_in->data.ushort_data, *r = vol_out->data.ushort_data;
    ulong width = vol_in->VolWidth, height = vol_in->VolHeight, depth = vol_in->VolDepth, x, z;
    
    for(z=0; z<depth; z++)
    {
       for(x=0; x<width; x++)
       {
          ErodeVerLine(f, width, height, k, g, h, h2, r);
          f++;
          r++;
       }
       f += width*(height-1);
       r += width*(height-1);
    }
} /* ImageGrayErodeVer */

void VolGrayDilateDepth(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out)
{
    ushort *f = vol_in->data.ushort_data, *r = vol_out->data.ushort_data;
    ulong width = vol_in->VolWidth, height = vol_in->VolHeight, depth = vol_in->VolDepth, x, y;

    for(y=0; y<height; y++)
    {
       for(x=0; x<width; x++)
       {
          DilateDepthLine(f, width, height, depth, k, g, h, h2, r);
          f++;
          r++;
       }
    }
} /* ImageGrayDilateDepth */


void VolGrayErodeDepth(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out)
{
    ushort *f = vol_in->data.ushort_data, *r = vol_out->data.ushort_data;
    ulong width = vol_in->VolWidth, height = vol_in->VolHeight, depth = vol_in->VolDepth, x, y;
    
    for(y=0; y<height; y++)
    {
       for(x=0; x<width; x++)
       {
          ErodeDepthLine(f, width, height, depth, k, g, h, h2, r);
          f++;
          r++;
       }
    }
} /* ImageGrayErodeDepth */
