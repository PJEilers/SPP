#ifndef _FILTER_H
#define _FILTER_H

#include "common.h"

/*typedef struct
{
    int self;
} ThreadFiltData;*/

void RefTreeAreaFilterBerger(double lambda, greyval_t *out, greyval_t *gval);
void Filter(GeneralThreadData *threadData, double lambda);
void RunFilter(GeneralThreadData *threadData, int nthreads);
void ParallelFilter(int nthreads, GeneralThreadData *threadData);
void *runfilt(void *arg);

#endif
