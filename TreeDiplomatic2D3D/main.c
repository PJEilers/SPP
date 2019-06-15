#include "common.h"
#include "quanttree.h"
#include "handleimages.h"
#include "radixsort.h"
#include "quantizedimage.h"
#include "refinetree.h"
#include "filter.h"

char *getExtension(const char *filename)
{
    char *dot = strrchr(filename, '.'); // A pointer to the last occurrence of character in str.
    return (char *)dot + 1;
}

void InitBuckets()
{
		NUMBUCKETS = 65536;
}

void printTimings()
{
	musec = (float)(end_sort - start_sort)/((float)tickspersec);
    printf("Sorting: %f s.\n", musec);
    
    musec = (float)(end_quimg - start_quimg)/((float)tickspersec);
    printf("Create Quantized Image: %f s.\n", musec);
    
    musec = (float)(end_qutree - start_qutree)/((float)tickspersec);
    printf("Quantized Tree: %f s.\n", musec);
    
    musec = (float)(end_ref - start_ref)/((float)tickspersec);
    printf("Refinement Tree: %f s.\n", musec);
    
    musec = (float)(end_filt - start_filt)/((float)tickspersec);
    printf("Filtering: %f s.\n", musec);
    
    musec = (float)(end_filt - start_sort)/((float)tickspersec);
    printf("Wall-Clock time: %f s.\n", musec);
}


int main (int argc, char *argv[])
{
    char *imgfname, *outfname = "out.fits";
    is3D = false;
    bitsPerPixel = 8; int r = 0;
    tickspersec = sysconf(_SC_CLK_TCK);  

    if (argc<7)
    {
		printf("Usage: %s <nthreads> <input image> <lambda> <bits per pixel> <output image> <is3D>.\n", argv[0]);
		exit(0);
    }
    
    nthreads = MIN(atoi(argv[1]), MAXTHREADS); // number of threads used for the parallel sorting algorithm and for the parallel building of the quantized max tree
    //inputQTZLEVELS = atoi(argv[2]); // number of quantized levels to work with
    int inputQTZLEVELS = nthreads; // atoi(argv[2]); // number of quantized levels to work with
    imgfname = argv[2];
    double lambda = atof(argv[3]);
    bitsPerPixel = atoi(argv[4]);
    outfname = argv[5];
    is3D = (atoi(argv[6]));
    width = height = depth = size = size2D = 0;
    ImageProperties img;	

    {if(!ReadFITS3D(imgfname, &img)){printf("fits image failed to read\n"); return(-1);}}

    size = width*height*depth;
    if(is3D)
      size2D = width*height;
    else
      size2D = width;
		
    
    GeneralThreadData *threadData = malloc(nthreads *sizeof(GeneralThreadData));
    
    int i;
    for(i = 0; i < nthreads; i++) {
       threadData[i].gval = gval; 
       threadData[i].nthreads = nthreads;
       threadData[i].img = img;
    }

    
    printf("Command: %s\n", argv[0]);
    printf("Filtering image '%s' using attribute area with lambda=%f.\n", imgfname, lambda);
    printf("FreeImage version: %s\n", FreeImage_GetVersion());
    printf("Image: Width=%ld Height=%ld Depth=%ld Size=%ld Size2D=%ld. ", width, height, depth, size, size2D);
    printf("nthreads for sorting and quantize the max tree: %d.\n", nthreads);
    printf("Size of int=%lu. Size of long=%ld. Size of float=%ld. Size of double=%ld.\n", sizeof(int), sizeof(long), sizeof(float), sizeof(double));
    
    //hmin = FindMin(gval); hmax = FindMax(gval);
    //printf("Min=%d. Max=%d.\n", hmin, hmax);
    
	/**************************************************************/ 
	printf("/*** Sort the pixels ***/\n");  
    InitBuckets();        
    start_sort = times(&tstruct);    
    // run radix sort
    RunRadixSortUnsigned(abs((int)ceil((double)bitsPerPixel/(16))));	
	end_sort = times(&tstruct);  
    
    /******************************************************************/     
    printf("/*** Calculate the quantized image ***/\n");
    start_quimg = times(&tstruct);
    CalculateQuantizedImage(inputQTZLEVELS, threadData);
	end_quimg = times(&tstruct);
    musec = (float)(end_quimg - start_quimg)/((float)tickspersec);
		
	/******************************************************************/ 	
        
    printf("/*** Build the max tree of the quantized image. (threads %d)***/\n", nthreads);    
    start_qutree = times(&tstruct);
    BuildMaxTreeOfQuantizedImage(threadData);
    end_qutree = times(&tstruct);   
    nthreadsRef = threadData->numQTZLEVELS;   
    printf("Pilot max-tree built. (numQTZLEVELS=%d).\n", threadData->numQTZLEVELS);
    
    /******************************************************************/ 	
    printf("/*** Refinement phase (threads %d) ***/\n", nthreadsRef);
	outRef =  malloc(size*sizeof(greyval_t));	
	RefineTreeBerger(nthreadsRef, threadData);
	printf("Refined max-tree built.\n"); 
	   
    /******************************************************************/
	printf("Init filtering\n");   
    start_filt = times(&tstruct);
    Filter(threadData, lambda);      
    end_filt = times(&tstruct);
    printf("End filtering\n");
    /******************************************************************/ 	
    printTimings();
	
	/**** Write output image   ****/   

    WriteFITS3D(outfname, imgfname, outRef);
    printf("Image written to '%s'\n", outfname);
	
    printf("\n");
    free(gval);
    free(SORTED);
	free(zpar);
	free(node_qu);
	free(node_ref);
	free(outRef);

    FreeImage_DeInitialise();
    return(0);
}
