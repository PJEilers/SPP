#ifndef _REFINETREE_H
#define _REFINETREE_H

#include "common.h"

/*typedef struct
{
    int self;
} ThreadRefData;*/

void RefineTreeBerger(int nthreads, GeneralThreadData *threadData);
void TreeAndCollectionsCreate(pixel_t size);
void RunRefinementThreads(GeneralThreadData *thdata, int nthreads);

// Manage parallel
//ThreadRefData *threfdata;
//ThreadRefData *MakeThreadRefData(int numthreads);
void FreeQuantized(GeneralThreadData *data, int numthreads);
void *rnc(void *arg);

//int GetNeighborsBerger(pixel_t p, pixel_t *neighbors);
int GetNeighborsBerger(pixel_t p, pixel_t *neighbors, pixel_t lwb, pixel_t upb);

pixel_t FINDROOT(pixel_t p);
#endif
