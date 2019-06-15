#include "common.h"
#include "sorting.h"
#include "quanttree.h"
#include "handleimages.h"

ThreadSortingData *MakeThreadSortingData(int numthreads)
{
    ThreadSortingData *data = malloc(numthreads *sizeof(ThreadSortingData));
    int i;

    for (i=0; i<numthreads; i++)
    {	
        data[i].self = i;
        data[i].hist = calloc(NUMBUCKETS, sizeof(pixel_t));
        //printf("Thread %d, LWB = %ld, UPB = %ld\n",i,LWB(i),UPB(i));
    }
    return(data);
}

void FreeThreadSortingData(ThreadSortingData *data, int numthreads)
{
    int i;
    for (i=0; i<numthreads; i++)
        free(data[i].hist);
    free(data);
}

void printHistogram(pixel_t *histogram)
{
	long i;
	printf("Histogram\n");
	for(i=0; i<NUMBUCKETS; i++)
	{
		if(histogram[i] != 0) printf("%ld) %ld\n", i, histogram[i]);
	}
	printf("---------\n");
}


void LocalHist(int self, greyval_t *gval, pixel_t *hist)
{
	pixel_t lwb, upb, i;
    lwb = LWB(self, nthreads);
    upb = UPB(self, nthreads);
    unsigned short radix;
    for (i=lwb; i<upb; ++i)
    {
		radix = *( (unsigned short*) &(gval[i]));
		hist[radix]++; // sort in ascending order
		
		//hist[gval[i]]++; // sort in ascending order
		//hist[ NUMBUCKETS -1 - gval[i] ]++; // sort in descending order
    }
    
    return;	
}

void CreateSortedArray(int self, greyval_t *gval, pixel_t *sorted, pixel_t *hist)
{
	pixel_t lwb, upb, i;
    lwb = LWB(self, nthreads);
    upb = UPB(self, nthreads);
    
    unsigned short radix;    
    i = lwb;
    for (i=lwb; i<upb; ++i)
    {
		//sorted[ hist[ gval[i] ]++ ] = gval[i]; // interested in px intensity rather than px coordinates
		//sorted[ hist[ NUMBUCKETS - 1 - gval[i] ]++ ] = i; // sort in descending order: remember to change also LocalHist function
		
		radix = *( (unsigned short*) &(gval[i]));
		sorted[ hist[ radix ]++ ] = i; // sort in ascending order
		//sorted[ hist[ gval[i] ]++ ] = i; // sort in ascending order
	}
    return;	
}

// Counting sort
void *csort(void *arg)
{	
    ThreadSortingData *thdata = (ThreadSortingData *) arg;
    int self = thdata->self;
    
    // build local histograms of each image partion
    LocalHist(self, gval, thdata->hist);    
    Barrier(self, nthreads);   
    
    // thread 0 collects the results from the other threads
    if(self == 0)
    {		
		int i;
		long j;
		histogram = calloc(NUMBUCKETS,  sizeof(pixel_t));
		
		// create the whole image histogram
		for(i=0; i<nthreads; i++)
		{
			ThreadSortingData *currThdata = (ThreadSortingData *) thsortingdata + i;
			for(j=0; j<NUMBUCKETS; j++)
			{
				histogram[j] += *(currThdata->hist+j);
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
			pixel_t sum = (thsortingdata)->hist[j], curr;
			for(i=1; i<nthreads; i++)
			{
				curr = (thsortingdata + i)->hist[j];
				(thsortingdata + i)->hist[j] = sum + histogram[j];
				sum += curr;
			}
		}
		((ThreadSortingData *) thsortingdata + 0)->hist = histogram;
		
		SORTED = (pixel_t *) malloc(width * depth * sizeof(pixel_t));
		
		//mapPixelToSorted = calloc(width * depth, sizeof(long));
	}
	
	// now that offset indexes are ready, every thread starts to create its part of the sorted array
	Barrier(self, nthreads);	
	
	
	// the threads now write the sorted array, each one starting from its initial offset
	CreateSortedArray(self, gval, SORTED, thdata->hist);	
	
    return NULL;
} 

void RunCountingSort(ThreadSortingData *allThdata, int nthreads)
{	
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, csort, (void *) (allThdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

void printSortedArray()
{
    pixel_t i;
    
    //pixel_t x,y;
    for(i=0; i<size; i++)
    {
        //x = SORTED[i]/width;
        //y = SORTED[i]%width;
        //printf("%ld) PixelPosition=(%ld, %ld)[%ld]; OrigGvalue=%d \n", i, y, x, SORTED[i], gval[SORTED[i]]);

        if( ((i>1) && (gval[SORTED[i]] < gval[SORTED[i-1]])) || ( (i>1) && (gval[SORTED[i]] == gval[SORTED[i-1]]) && (SORTED[i] < SORTED[i-1]) ) )
        {
            printf("ERROR\n"); // in the IF, use '>' if descending ordering is used; '<' otherwise
            exit(0);
        }
    }
    
    printf("===========================================================\n");
}

void ComputeHistogram()
{
    long i;
    unsigned short radix;
    pixel_t myhist[NUMBUCKETS];
    for(i=0; i<NUMBUCKETS; i++)
    {
        myhist[i] = 0;
    }

    for(i=0; i<size; i++)
    {
		radix = *( (unsigned short*) &(gval[i]));
		myhist[radix]++;
        //myhist[(gval[i])]++;
    }

    printHistogram(myhist);
    return;
}

void SortImagePixels()
{
	printf("Counting sort\n");
	
    // Sort the pixels in array SORTED
	thsortingdata = MakeThreadSortingData(nthreads);
	RunCountingSort(thsortingdata, nthreads);
	FreeThreadSortingData(thsortingdata, nthreads);	

	//printSortedArray();
	//ComputeHistogram();
}

