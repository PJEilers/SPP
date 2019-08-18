#ifndef _TREEDIPLOMATIC_H
#define _TREEDIPLOMATIC_H

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
#define NUMBUCKETS 65536
#define CONNECTIVITY 6 // keep it at 6 so it is ready for both 2D and 3D. Code realised if it dealing with a 2D or 3D image.

/*********************************************/
#define USEFLOATPOINT 1
typedef float greyval_t; //intensities
/*********************************************/


typedef long pixel_t; //coordinates
#define USEFLOATPOINT_ONLYPOSITIVE 0
typedef short bool;
typedef unsigned char ubyte;


#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))
#define LWB(self, nthreads, size2D, depth) (size2D*(((self)*depth)/nthreads))
#define UPB(self, nthreads, size2D, depth) (size2D*(((self+1)*depth)/nthreads))
#define bottom (-1)

//int nthreadsRef; // number of threads of the NajCou refinement, equals to the number of quantized grey levels numQTZLEVELS


typedef struct ImageProperties 
{
    long width, height, depth, size, size2D;  /* precondition: width <= size/nthreads */
    short bitsPerPixel;
} ImageProperties;

typedef struct Queue
{
    pixel_t *Pixels;
    pixel_t Head;
    pixel_t Tail; /* First free place in queue, empty if Head=Tail */
} Queue;

typedef struct MaxNode
{
    pixel_t parent;
    long Area;
} MaxNode;


typedef struct Node Node;
struct Node{
	greyval_t filter;
	bool isFiltered;
    pixel_t parent;
    long Area;
};

typedef struct GeneralThreadData
{
    int self;
    int numQTZLEVELS;
    double lambda;
    Queue *thisqueue;
    greyval_t *gval;
    greyval_t *outRef;
    int *gval_qu; 
    int nthreads;
    pixel_t *pxStartPosition;
    pixel_t *pxEndPosition;
    ImageProperties img;
    pixel_t *SORTED;
    bool *reached_qu;
    MaxNode *node_qu;
    Node *node_ref;
    pixel_t *zpar;
} GeneralThreadData;

/***** Timings ****/

typedef struct Timings
{
    clock_t start;
    struct tms tstruct;
    float musec;
    long tickspersec;

    clock_t start_sort, end_sort; 
    clock_t start_quimg, end_quimg;
    clock_t start_qutree, end_qutree;
    clock_t start_ref, end_ref;
    clock_t start_filt, end_filt;
} Timings;

void FreeThreadData(GeneralThreadData *threadData);

/*Quantize Image Functions*/

int CreateQuantizedImage(int inputQTZLEVELS, GeneralThreadData *threadData);
void CalculateQuantizedImage(int inputQTZLEVELS, GeneralThreadData *threadData);
int WorkOutBoundaries(int inputQTZLEVELS, GeneralThreadData *threadData);

void CreateQuantizedImageInParallel(int numQtzLevsUsed, GeneralThreadData *threadData);
void MakeThreadPixelPositions(int numthreads, GeneralThreadData *data);
void MakeThreadQntImgData(int numthreads, GeneralThreadData *threadData);
void RunThreadsWritingQuantizedImage(GeneralThreadData *thdata, int nthreads);
void *wqi(void *arg);

/* Quanttree Functions */

#define QueueFirst(hq,h)  (hq[h].Pixels[hq[h].Head++])
#define QueueAdd(hq,h,p)  hq[h].Pixels[hq[h].Tail++] = p
#define QueueNotEmpty(hq,h)  (hq[h].Head != hq[h].Tail)



void *SafeMalloc(int n);
void *SafeCalloc(int nmemb, pixel_t size);
void Psa(int p);
void Vsa(int p);
void Barrier(int self, int nthreads);
pixel_t levroot(pixel_t x, int *gvalues, MaxNode *node_qu);
pixel_t Par(pixel_t x, int *gvalues, MaxNode *node_qu);
void levrootfix(pixel_t lwb, pixel_t upb, int *gvalues, MaxNode *node_qu);
 
Queue *QueueCreate(pixel_t imgsize, int numQTZLEVELS);
void FillPixelsInQueue(Queue *hq, pixel_t *numpixelsperlevel, int numQTZLEVELS);
void QueueDelete(Queue *hq);

int GetNeighbors(pixel_t p, pixel_t x, pixel_t y, pixel_t z, pixel_t *neighbors, pixel_t lwb, pixel_t upb, ImageProperties img);
int LocalTreeFlood(int self,  Queue *set, pixel_t *lero, int lev, long *thisarea, MaxNode *nodes, int *gvalues, GeneralThreadData *threadData);
void Connect(pixel_t x, pixel_t y, MaxNode *node, int *gvalues, greyval_t *gval);
void Fuse(int self, int i, MaxNode *node, int *gvalues, ImageProperties img, greyval_t *gval, int nthreads);
void MakeThreadData(int numthreads, GeneralThreadData *data);
void FreeThreadQueues(GeneralThreadData *data, int nthreads);
void BuildQuantizedTree(GeneralThreadData *thdata, int nthreads);
void *ccaf(void *arg);
 
int BuildMaxTreeOfQuantizedImage();


/* Radix Sort Functions*/

typedef struct ThreadSortingDataRS ThreadSortingDataRS;

struct ThreadSortingDataRS 
{
    int self;
    pixel_t *hist;
    int numStepsRS;
    ImageProperties img;
    pixel_t *SORTEDRS[2];
    pixel_t *HISTOGRAMRS[2];
    greyval_t *gval;
    int nthreads;
    ThreadSortingDataRS *thsortingdataRS;
};

pixel_t *RunRadixSortUnsigned(int numSteps, ImageProperties img, greyval_t *gval, int nthreads);
void RunRSCountingSort(ThreadSortingDataRS *allThdata, int nthreads);
//void FreeThreadSortingDataRS (ThreadSortingDataRS *data, int nthreads);
void *csortRS(void *arg);
void *fpx(void *arg);
void *fpxback(void *arg);
void RunThreadsFlipping(ThreadSortingDataRS *thdata, int nthreads);
void RunThreadsFlippingBack(ThreadSortingDataRS *thdata, int nthreads);
void LocalHistRS(int self, pixel_t *SORTED, pixel_t *hist, int step, ThreadSortingDataRS thdata);
ThreadSortingDataRS *MakeThreadRadixSortData(int numthreads, int numSteps, ImageProperties img, greyval_t *gval);
void FreeSortingDataRS(ThreadSortingDataRS *data, int numthreads, int numSteps);
void CreateSortedArrayRS(int self, pixel_t *sortedNew, pixel_t *sortedOld, pixel_t *hist, int step, ThreadSortingDataRS *thdata);

void RefineTreeBerger(int nthreads, GeneralThreadData *threadData);
void TreeAndCollectionsCreate(pixel_t size, GeneralThreadData *threadData);
void RunRefinementThreads(GeneralThreadData *thdata, int nthreads);

/* Refine Tree Functions */

// Manage parallel
//ThreadRefData *MakeThreadRefData(int numthreads);
void FreeQuantized(GeneralThreadData *data, int numthreads);
void *rnc(void *arg);

//int GetNeighborsBerger(pixel_t p, pixel_t *neighbors);
int GetNeighborsBerger(pixel_t p, pixel_t *neighbors, pixel_t lwb, pixel_t upb, ImageProperties img);
pixel_t DescendRoots(pixel_t q, int myLev, int *gval_qu, MaxNode *node_qu);

pixel_t FINDROOT(pixel_t p, pixel_t *zpar);

/* Filtering Functions */

void RefTreeAreaFilterBerger(double lambda, greyval_t *out, greyval_t *gval, long size, Node *node_ref);
void Filter(GeneralThreadData *threadData, double lambda);
void RunFilter(GeneralThreadData *threadData, int nthreads);
void ParallelFilter(int nthreads, GeneralThreadData *threadData);
void *runfilt(void *arg);


#endif
