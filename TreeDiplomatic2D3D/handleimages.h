#ifndef _HANDLEIMAGES_H
#define _HANDLEIMAGES_H

#include <FreeImage.h>
#include <fitsio.h>
#include <math.h>
#include "TreeDiplomatic.h"

FIBITMAP* GenericLoader(const char* lpszPathName, int flag);
greyval_t *ReadFITS(char *filename, ImageProperties *img);
greyval_t FindMin(greyval_t *gvalues);
greyval_t FindMax(greyval_t *gval);
short ReadTIFF(char *fname);
int WriteFITS(char *filename, char *inputimagefilename, greyval_t *out);
int WriteTIFF(char *fname, greyval_t *img);
void WriteTIFF_2(char *fname, greyval_t *img, pixel_t width, pixel_t height, int bitspp);
greyval_t *ReadFITS3D(char *filename, ImageProperties *img);
int WriteFITS3D(char *filename, char *inputimagefilename, greyval_t *out, ImageProperties img);
int ReadFITS3D_CUBE(char *filename);
int WriteFITS3DCUBE(char *filename, char *inputimagefilename, greyval_t *out);

#endif
