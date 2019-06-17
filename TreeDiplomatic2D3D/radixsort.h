#ifndef _RADIXSORT_H
#define _RADIXSORT_H

typedef struct
{
    int self;
    pixel_t *hist;
    int numStepsRS;
    ImageProperties img;
    pixel_t *SORTEDRS[2];
    pixel_t *HISTOGRAMRS[2];
    greyval_t *gval;
} ThreadSortingDataRS;

typedef struct
{
    int self;
    greyval_t *gval;
} ThreadFlipData;


ThreadSortingDataRS *thsortingdataRS;
pixel_t *HISTOGRAMRS[2];

pixel_t *RunRadixSortUnsigned(int numSteps, ImageProperties img, greyval_t *gval);
void RunRSCountingSort(ThreadSortingDataRS *allThdata, int nthreads);
void *csortRS(void *arg);
void *fpx(void *arg);
void *fpxback(void *arg);
void RunThreadsFlipping(ThreadSortingDataRS *thdata, int nthreads);
void RunThreadsFlippingBack(ThreadSortingDataRS *thdata, int nthreads);
void LocalHistRS(int self, pixel_t *SORTED, pixel_t *hist, int step, ThreadSortingDataRS thdata);
ThreadSortingDataRS *MakeThreadRadixSortData(int numthreads, int numSteps, ImageProperties img, greyval_t *gval);
void FreeThreadSortingDataRS(ThreadSortingDataRS *data, int numthreads);
void CreateSortedArrayRS(int self, pixel_t *sortedNew, pixel_t *sortedOld, pixel_t *hist, int step, ThreadSortingDataRS *thdata);

#endif
