#include "common.h"
#include "quantizedimage.h"
#include "quanttree.h"


void MakeThreadPixelPositions(int numthreads, GeneralThreadData *data) {
    
    pixel_t *pxStartPosition = calloc(numthreads, sizeof(pixel_t));
    pixel_t *pxEndPosition = calloc(numthreads, sizeof(pixel_t));
    int i;
    for (i=0; i<numthreads; i++)
    {
        data[i].pxStartPosition = pxStartPosition;
        data[i].pxEndPosition = pxEndPosition;
        data[i].numQTZLEVELS = numthreads;
    }
}

void MakeThreadQntImgData(int numthreads, GeneralThreadData *data)
{
    /* create the quantized image */
    long size = data->img.size;
    int *gval_qu = malloc(size*sizeof(int));

    if (gval_qu==NULL)
    {
        fprintf(stderr, "out of memory! \n");
        free(gval_qu);
        exit(0);
    }
    int i;

    for (i=0; i<numthreads; i++)
    {
        data[i].gval_qu = gval_qu;
        data[i].self=i;
        data[i].numQTZLEVELS = numthreads;
    }
}


void RunThreadsWritingQuantizedImage(GeneralThreadData *thdata, int nthreads)
{
    pthread_t threadID[MAXTHREADS];
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, wqi, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

/* CreateQuantizedImage */
int CreateQuantizedImage(int inputQTZLEVELS, GeneralThreadData *threadData)
{	
    long size = threadData->img.size;
    pixel_t numRemPixels = size;
    pixel_t expectedNumPixelPerQtzLev = numRemPixels/(inputQTZLEVELS);
    pixel_t *pxStartPosition = threadData->pxStartPosition;
    pixel_t *pxEndPosition = threadData-> pxEndPosition;
    pixel_t *SORTED = threadData->SORTED;
    int *gval_qu = threadData->gval_qu;
    greyval_t *gval = threadData->gval;
    int currentQtzLev = 0;
    pixel_t px=0, count = 0;
    greyval_t prevBoundGval = 0;
    int idxThread = 0;

    while(px<size)
    {
		pxStartPosition[idxThread]	= px;
		
		// iterate until the optimal number of pixels for the current thread
        while(px<size && (count < expectedNumPixelPerQtzLev || currentQtzLev == inputQTZLEVELS-1))
        {
            gval_qu[SORTED[px]] = currentQtzLev;
            count++;
            prevBoundGval = gval[SORTED[px]];
            px++;
        }
        
        // iterate until the next intensity change to map different intensities on different threads    
        while(px<size && prevBoundGval == gval[SORTED[px]])
        {
            gval_qu[SORTED[px]] = currentQtzLev;
            count++;
            px++;
        }
        
        // update thread info
        pxEndPosition[idxThread] = px-1;
                
        // work out the optimal number of pixels for the next thread
        numRemPixels = numRemPixels - count;
        expectedNumPixelPerQtzLev = numRemPixels/(inputQTZLEVELS-currentQtzLev);
        
        currentQtzLev++;        
        idxThread++;
        count=0;
     }

    return currentQtzLev; // return the number of quantization levels used
}  

int WorkOutBoundaries(int inputQTZLEVELS, GeneralThreadData *threadData)
{
    long size = threadData->img.size;
    pixel_t numRemPixels = size;
    pixel_t expectedNumPixelPerQtzLev = numRemPixels/(inputQTZLEVELS);
    pixel_t *pxStartPosition = threadData->pxStartPosition;
    pixel_t *pxEndPosition = threadData-> pxEndPosition;
    pixel_t *SORTED = threadData->SORTED;
    
    greyval_t *gval = threadData->gval;
    //printf("Optimal number of pixels per thread = %ld\n", expectedNumPixelPerQtzLev);

    int currentQtzLev = 0;
    pixel_t px=0;
    greyval_t prevBoundGval = 0;
    int idxThread = 0;

    while(px<size)
    {
		pxStartPosition[idxThread]	= px;
		
		px += (expectedNumPixelPerQtzLev - 1); // point to the last gvalue of that thread.
		
		if(px<size)
		{
			prevBoundGval = gval[SORTED[px]];
			px++;
		
			while( (px<size) && ( (prevBoundGval == gval[SORTED[px]]) || (currentQtzLev == inputQTZLEVELS-1) ) )
			{			
				px++;
			}
		
			pxEndPosition[idxThread] = px - 1;
		}
	    else 
			pxEndPosition[idxThread] = size-1;
        
        currentQtzLev++;        
        idxThread++;
     }
     
    return currentQtzLev; // return the number of quantization levels used
}

void *wqi(void *arg)
{	
    GeneralThreadData *threfdata = (GeneralThreadData *) arg;
    int self = threfdata->self;
    pixel_t *pxStartPosition = threfdata->pxStartPosition;
    pixel_t *pxEndPosition = threfdata->pxEndPosition;
    pixel_t *SORTED = threfdata->SORTED;
    int *gval_qu = threfdata->gval_qu;
    
    pixel_t i;
    for(i=pxStartPosition[self]; i <= pxEndPosition[self]; i++)
    {
		gval_qu[SORTED[i]] = self;
	}	   
    return NULL;
}

void CreateQuantizedImageInParallel(int numQtzLevsUsed, GeneralThreadData *threadData)
{
    MakeThreadQntImgData(numQtzLevsUsed, threadData);    
    RunThreadsWritingQuantizedImage(threadData, numQtzLevsUsed);
}

void CalculateQuantizedImage(int inputQTZLEVELS, GeneralThreadData *threadData)
{
    


    
        MakeThreadPixelPositions(inputQTZLEVELS, threadData);
        
    
		// First Way: in parallel
		int numQTZLEVELS = WorkOutBoundaries(inputQTZLEVELS, threadData);

		if(numQTZLEVELS != inputQTZLEVELS)
		{
			printf("Trying second method and sacrifice load balance\n");
			// Second Way: use the old sequential (worst load balance)
			numQTZLEVELS = CreateQuantizedImage(inputQTZLEVELS, threadData);
			if(numQTZLEVELS != inputQTZLEVELS)
			{
				printf("ERROR: input quantized levels = %d, but actually used %d ! ", inputQTZLEVELS, numQTZLEVELS);
				printf("Exiting: please reduce the number of threads!\n");
				//printf("Sacrifice load balance and maintain the input number of quantization levels using e.g. \"CreateQuantizedImage()\".\n");
				exit(0);
			}
			return;
		}		
		// if everything goes good:
		CreateQuantizedImageInParallel(numQTZLEVELS, threadData);	
    return;
}
