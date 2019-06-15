#ifndef _SORTING_H
#define _SORTING_H

typedef struct
{
    int self;
    pixel_t *hist;
} ThreadSortingData;

ThreadSortingData *thsortingdata;

// Used for parallel counting sort
pixel_t *histogram;

ThreadSortingData *MakeThreadSortingData(int numthreads);
void FreeThreadSortingData(ThreadSortingData *data, int numthreads);
void printHistogram(pixel_t *histogram);
void LocalHist(int self, greyval_t *gval, pixel_t *hist);
void CreateSortedArray(int self, greyval_t *gval, pixel_t *sorted, pixel_t *hist);
void *csort(void *arg);
void RunCountingSort(ThreadSortingData *allThdata, int nthreads);
void printSortedArray();
void ComputeHistogram();
void SortImagePixels();
#endif
