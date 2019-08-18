#include "TreeDiplomatic.h"

void RunRSCountingSort(ThreadSortingDataRS *allThdata, int nthreads)
{	
    pthread_t threadID[MAXTHREADS];
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, csortRS, (void *) (allThdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

void FreeSortingDataRS(ThreadSortingDataRS *data, int numthreads, int numSteps) {
    int i;
    free(data->SORTEDRS[(numSteps+1)%2]);
    free(data->HISTOGRAMRS[(numSteps)%2]);
    for (i=0; i<numthreads; i++)
        free(data[i].hist);
    free(data);
}

ThreadSortingDataRS *MakeThreadRadixSortData(int numthreads, int numSteps, ImageProperties img, greyval_t *gval)
{
    ThreadSortingDataRS *data = malloc(numthreads *sizeof(ThreadSortingDataRS));
    pixel_t *SORTEDRS[2];
    pixel_t *HISTOGRAMRS[2];
    HISTOGRAMRS[0] = (pixel_t *) malloc(NUMBUCKETS*sizeof(pixel_t));
    HISTOGRAMRS[1] = (pixel_t *) malloc(NUMBUCKETS*sizeof(pixel_t));
    SORTEDRS[0] = (pixel_t *) malloc(img.size * sizeof(pixel_t));
    SORTEDRS[1] = (pixel_t *) malloc(img.size * sizeof(pixel_t));
    int i;
    
    for (i=0; i<numthreads; i++)
    {	
        data[i].self = i;
        data[i].hist = calloc(NUMBUCKETS, sizeof(pixel_t));	
        data[i].numStepsRS = numSteps;
        data[i].img = img;
        data[i].gval = gval;
        data[i].nthreads = numthreads;
        data[i].HISTOGRAMRS[0] = HISTOGRAMRS[0];
        data[i].HISTOGRAMRS[1] = HISTOGRAMRS[1];
        data[i].SORTEDRS[0] = SORTEDRS[0];
        data[i].SORTEDRS[1] = SORTEDRS[1];

        //printf("Thread %d, LWB = %ld, UPB = %ld\n",i,LWB(i, numthreads),UPB(i, numthreads));
    }
    data[0].thsortingdataRS = data;
    return(data);
}

/* works on 2 bytes rather than 1 byte each loop */
void LocalHistRS2(int self, pixel_t *sorted, pixel_t *hist, int step, ThreadSortingDataRS *thdata)
{
	pixel_t lwb, upb, i;
    int nthreads = thdata->nthreads;
    lwb = LWB(self, nthreads, thdata->img.size2D, thdata->img.depth);
    upb = UPB(self, nthreads, thdata->img.size2D, thdata->img.depth);
    unsigned short radix;
    greyval_t *gval = thdata->gval;
        
    memset(hist, 0, sizeof(pixel_t)*NUMBUCKETS);
    
    if(step == 0)
    {
		for (i=lwb; i<upb; ++i)
		{
			radix = *( (unsigned short*) &(gval[i]));			
			hist[radix]++;
		}
	}
	else
	{
		for (i=lwb; i<upb; ++i)
		{			
			radix = *( (unsigned short*) (&(gval[sorted[i]]))+step);			
			hist[radix]++;
		}
	}
	    
    return;	
}

void CreateSortedArrayRS2(int self, pixel_t *sortedNew, pixel_t *sortedOld, pixel_t *hist, int step, ThreadSortingDataRS *thdata)
{
	pixel_t lwb, upb, i;
    int nthreads = thdata->nthreads;
    lwb = LWB(self, nthreads, thdata->img.size2D, thdata->img.depth);
    upb = UPB(self, nthreads, thdata->img.size2D, thdata->img.depth);
    greyval_t *gval = thdata->gval;
    unsigned short radix;
    
    if(step == 0)
    {
		for (i=lwb; i<upb; ++i)
		{	
			radix = *( (unsigned short*) &(gval[i]));
			sortedNew[ hist[ radix ]++ ] = i; // sort in ascending order			
		}
	}
	else
	{
		for (i=lwb; i<upb; ++i)
		{		
			radix = *( (unsigned short*) (&(gval[sortedOld[i]]))+step);
			sortedNew[ hist[ radix ]++ ] = sortedOld[i]; // sort in ascending order
		}		
	}
    return;	
}

// counting sort for radix sort
void *csortRS(void *arg)
{	
    ThreadSortingDataRS *thdata = (ThreadSortingDataRS *) arg;
    int self = thdata->self;
    int numStepsRS = thdata->numStepsRS;
    int nthreads = thdata->nthreads;
    pixel_t *histogram;
    pixel_t *HISTOGRAMRS[2];
    pixel_t *SORTEDRS[2];
    HISTOGRAMRS[0] = thdata->HISTOGRAMRS[0];
    HISTOGRAMRS[1] = thdata->HISTOGRAMRS[1];
    SORTEDRS[0] = thdata->SORTEDRS[0];
    SORTEDRS[1] = thdata->SORTEDRS[1];
    int step;
    for(step=0; step<numStepsRS; step++)
    {  
		// build local histograms of each image partion
		LocalHistRS2(self, SORTEDRS[step%2], thdata->hist, step, thdata);
		Barrier(self, nthreads);

		// thread 0 collects the results from the other threads
		if(self == 0)
		{
                        ThreadSortingDataRS *thsortingdataRS = thdata->thsortingdataRS;
			pixel_t i,j;
			histogram = HISTOGRAMRS[step%2];
			
			memset(histogram, 0, sizeof(pixel_t)*NUMBUCKETS);
			
			// create the whole image histogram
			for(i=0; i<nthreads; i++)
			{
				ThreadSortingDataRS *currThdataRS = (ThreadSortingDataRS *) thsortingdataRS + i;
				for(j=0; j<NUMBUCKETS; j++)
				{
					histogram[j] += *(currThdataRS->hist+j);
				}
			}
			
			// overwrite 'histogram' with the initial offset indexes of each intensity value
			pixel_t prevhist, total= 0;
			for(j=0; j<NUMBUCKETS; j++)
			{
				prevhist = histogram[j];				
				histogram[j] = total;
				total += prevhist;
			}
			
			// overwrite the 'hist' structure in every thread and write the proper initial offset
			for(j=0; j<NUMBUCKETS; j++)
			{
				pixel_t sum = (thsortingdataRS)->hist[j], curr;
				for(i=1; i<nthreads; i++)
				{
					curr = (thsortingdataRS + i)->hist[j];
					(thsortingdataRS + i)->hist[j] = sum + histogram[j];
					sum += curr;
				}
			}
			
			((ThreadSortingDataRS *) thsortingdataRS + 0)->hist = histogram;
			//memcpy(((ThreadSortingDataRS *) thsortingdataRS + 0)->hist, histogram, NUMBUCKETS*sizeof(pixel_t) );
		}
		
		// now that offset indexes are ready, every thread starts to create its part of the sorted array
		Barrier(self, nthreads);
		CreateSortedArrayRS2(self, SORTEDRS[((step+1)%2)], SORTEDRS[step%2], thdata->hist, step, thdata);
		Barrier(self, nthreads);
	}
	
	return NULL;
}

void FloatFlip(unsigned int *f)
{
	unsigned int mask = - (unsigned int)(*f >> 31) | 0x80000000;
	*f ^= mask;
	return;
}

void IFloatFlip(unsigned int *f)
{
	unsigned int mask = ((*f >> 31) - 1) | 0x80000000;
	*f ^= mask;
	return; // *f ^ mask;
}

void DoubleFlip(unsigned long *d)
{
	unsigned long mask = - (unsigned long)(*d >> 63) | 0x8000000000000000;
	*d ^= mask;
	return;
}

void IDoubleFlip(unsigned long *d)
{
	unsigned long mask = ((*d >> 63) - 1) | 0x8000000000000000;
	*d ^= mask;
	return;// *d ^ mask;
}

//flip pixels
void *fpx(void *arg)
{
    ThreadSortingDataRS *threfdata = (ThreadSortingDataRS *) arg;
    greyval_t *gval = threfdata->gval;
    int self = threfdata->self;
    short bitsPerPixel = threfdata->img.bitsPerPixel;
    int nthreads = threfdata-> nthreads;
    pixel_t lwb, upb, i;
    lwb = LWB(self, nthreads, threfdata->img.size2D, threfdata->img.depth);
    upb = UPB(self, nthreads, threfdata->img.size2D, threfdata->img.depth);
	
	if(bitsPerPixel == 32)
	{
		for(i=lwb; i<upb; i++)
		{
			FloatFlip((unsigned int *)&gval[i]);
		}
	}
	else if(bitsPerPixel == 64)
	{
		for(i=lwb; i<upb; i++)
		{
			DoubleFlip((unsigned long *)&gval[i]);			
		}		
	}
	
	//Barrier(self, nthreads);		   
    return NULL;
}

//flip pixels back
void *fpxback(void *arg)
{
    ThreadSortingDataRS *threfdata = (ThreadSortingDataRS *) arg;
    greyval_t *gval = threfdata->gval;
    short bitsPerPixel = threfdata->img.bitsPerPixel;
    int self = threfdata->self;
    int nthreads = threfdata-> nthreads;
    pixel_t lwb, upb, i;
    lwb = LWB(self, nthreads, threfdata->img.size2D, threfdata->img.depth);
    upb = UPB(self, nthreads, threfdata->img.size2D, threfdata->img.depth);
	
	if(bitsPerPixel == 32)
	{
		for(i=lwb; i<upb; i++)
		{
			IFloatFlip((unsigned int *)&gval[i]);
		}
	}
	else if(bitsPerPixel == 64)
	{
		for(i=lwb; i<upb; i++)
		{
			IDoubleFlip((unsigned long *)&gval[i]);			
		}		
	}
	
	//Barrier(self, nthreads);		   
    return NULL;
}

// ThreadFlipData *MakeThreadFlipData(int numthreads)
// {
//     ThreadFlipData *data = malloc(numthreads *sizeof(ThreadFlipData));
//     int i;
// 
//     for (i=0; i<numthreads; i++)
//     {
//         data[i].self=i;
//     }
//     return(data);
// }
// 
// void FreeThreadFlipData(ThreadFlipData *data, int numthreads)
// {
//     free(data);
// }

void RunThreadsFlipping(ThreadSortingDataRS *thdata, int nthreads)
{
    pthread_t threadID[MAXTHREADS];
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, fpx, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

void RunThreadsFlippingBack(ThreadSortingDataRS *thdata, int nthreads)
{
    pthread_t threadID[MAXTHREADS];
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, fpxback, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

/*void FlipParallel(int nthreads, greyval_t *gval)
{
    ThreadFlipData *thflipdata = MakeThreadFlipData(nthreads);    
    RunThreadsFlipping(thflipdata, nthreads);
    FreeThreadFlipData(thflipdata, nthreads);
}

void FlipBackParallel(int nthreads, greyval_t *gval)
{
    ThreadFlipData *thflipdata = MakeThreadFlipData(nthreads, gval);    
    RunThreadsFlippingBack(thflipdata, nthreads);
    FreeThreadFlipData(thflipdata, nthreads);
}*/

pixel_t *RunRadixSortUnsigned(int numSteps, ImageProperties img, greyval_t *gval, int nthreads)
{	
    
	printf("Radix Sort (steps=%d)\n", numSteps);
	
	ThreadSortingDataRS *thsortingdataRS = MakeThreadRadixSortData(nthreads, numSteps, img, gval);
	
	
	if(USEFLOATPOINT && !USEFLOATPOINT_ONLYPOSITIVE)
	{
                RunThreadsFlipping(thsortingdataRS, nthreads);
	}
	
	RunRSCountingSort(thsortingdataRS, nthreads);	
	pixel_t *SORTED = thsortingdataRS->SORTEDRS[numSteps%2];

	if(USEFLOATPOINT && !USEFLOATPOINT_ONLYPOSITIVE)
	{
                RunThreadsFlippingBack(thsortingdataRS,
                                       nthreads);
	}
	
	FreeSortingDataRS(thsortingdataRS, nthreads, numSteps);
	return SORTED;	
}


