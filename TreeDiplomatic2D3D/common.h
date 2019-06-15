#ifndef _COMMON_H
#define _COMMON_H

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <fitsio.h>
#include <math.h>
#include <FreeImage.h>


#define false 0
#define true  1
#define MAXTHREADS 128
#define CONNECTIVITY 6 // keep it at 6 so it is ready for both 2D and 3D. Code realised if it dealing with a 2D or 3D image.

/*********************************************/
#define USEFLOATPOINT 1
typedef float greyval_t; //intensities
/*********************************************/


typedef long pixel_t; //coordinates
#define USEFLOATPOINT_ONLYPOSITIVE 0
typedef short bool;
bool is3D;

greyval_t hmin; 
greyval_t hmax;

typedef unsigned char ubyte;


#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))
#define LWB(self, nthreads) (size2D*(((self)*depth)/nthreads))
#define UPB(self, nthreads) (size2D*(((self+1)*depth)/nthreads))
#define bottom (-1)

int nthreads;
int nthreadsRef; // number of threads of the NajCou refinement, equals to the number of quantized grey levels numQTZLEVELS

long width, height, depth, size, size2D;  /* precondition: width <= size/nthreads */

greyval_t *gval;
greyval_t *outRef;

short bitsPerPixel;

typedef struct ImageProperties 
{
    long width, height, depth, size, size2D;  /* precondition: width <= size/nthreads */
} ImageProperties;

typedef struct Queue
{
    pixel_t *Pixels;
    pixel_t Head;
    pixel_t Tail; /* First free place in queue, empty if Head=Tail */
} Queue;

typedef struct GeneralThreadData
{
    int self;
    int numQTZLEVELS;
    double lambda;
    Queue *thisqueue;
    greyval_t *gval;
    int *gval_qu; 
    int nthreads;
    pixel_t *pxStartPosition;
    pixel_t *pxEndPosition;
    ImageProperties img;
} GeneralThreadData;

typedef struct MaxNode
{
    pixel_t parent;
    long Area;
} MaxNode;

MaxNode *node_qu;
bool *reached_qu;

// Used for parallel counting sort
pixel_t *SORTED;
pixel_t *SORTEDRS[2];
unsigned int NUMBUCKETS;



typedef struct Node Node;
struct Node{
	greyval_t filter;
	bool isFiltered;
    pixel_t parent;
    long Area;
};

Node *node_ref;
pixel_t *zpar;

/***** Timings ****/
clock_t start;
struct tms tstruct;
float musec;
long tickspersec;

clock_t start_sort, end_sort; 
clock_t start_quimg, end_quimg;
clock_t start_qutree, end_qutree;
clock_t start_ref, end_ref;
clock_t start_filt, end_filt;

#endif