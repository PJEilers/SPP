#include "TreeDiplomatic.h"
#include "handleimages.h"

char *getExtension(const char *filename)
{
    char *dot = strrchr(filename, '.'); // A pointer to the last occurrence of character in str.
    return (char *)dot + 1;
}

void FreeThreadData(GeneralThreadData *threadData) {
    free(threadData->gval);
    free(threadData->outRef);
    free(threadData->gval_qu);
    free(threadData->pxStartPosition);
    free(threadData->pxEndPosition);
    free(threadData->SORTED);
    free(threadData->reached_qu);
    free(threadData->node_qu);
    free(threadData->node_ref);
    free(threadData->zpar);
    free(threadData);
}

void printTimings(Timings timings)
{
    timings.musec = (float)(timings.end_sort - timings.start_sort)/((float)timings.tickspersec);
    printf("Sorting: %f s.\n", timings.musec);
    
    timings.musec = (float)(timings.end_quimg - timings.start_quimg)/((float)timings.tickspersec);
    printf("Create Quantized Image: %f s.\n", timings.musec);
    
    timings.musec = (float)(timings.end_qutree - timings.start_qutree)/((float)timings.tickspersec);
    printf("Quantized Tree: %f s.\n", timings.musec);
    
    timings.musec = (float)(timings.end_ref - timings.start_ref)/((float)timings.tickspersec);
    printf("Refinement Tree: %f s.\n", timings.musec);
    
    timings.musec = (float)(timings.end_filt - timings.start_filt)/((float)timings.tickspersec);
    printf("Filtering: %f s.\n", timings.musec);
    
    timings.musec = (float)(timings.end_filt - timings.start_sort)/((float)timings.tickspersec);
    printf("Wall-Clock time: %f s.\n", timings.musec);
}


int main (int argc, char *argv[])
{
    char *imgfname, *outfname = "out.fits";
    bool is3D = false;
    short bitsPerPixel = 8; 
    Timings timings;
    timings.tickspersec = sysconf(_SC_CLK_TCK); 

    if (argc<7)
    {
		printf("Usage: %s <nthreads> <input image> <lambda> <bits per pixel> <output image> <is3D>.\n", argv[0]);
		exit(0);
    }
    
    int nthreads = MIN(atoi(argv[1]), MAXTHREADS); // number of threads used for the parallel sorting algorithm and for the parallel building of the quantized max tree
    //inputQTZLEVELS = atoi(argv[2]); // number of quantized levels to work with
    int inputQTZLEVELS = nthreads; // atoi(argv[2]); // number of quantized levels to work with
    imgfname = argv[2];
    double lambda = atof(argv[3]);
    bitsPerPixel = atoi(argv[4]);
    outfname = argv[5];
    is3D = (atoi(argv[6]));
    //width = height = depth = size = size2D = 0;
    ImageProperties img;	
    greyval_t *gval = ReadFITS3D(imgfname, &img);
   // {if(!ReadFITS3D(imgfname, &img){printf("fits image failed to read\n"); return(-1);}}
    img.size = img.width*img.height*img.depth;
    if(is3D)
      img.size2D = img.width*img.height;
    else
      img.size2D = img.width;
		
    img.bitsPerPixel = bitsPerPixel;
    GeneralThreadData *threadData = malloc(nthreads *sizeof(GeneralThreadData));
    
    
    int i;
    greyval_t *outRef =  malloc(img.size*sizeof(greyval_t));	
    for(i = 0; i < nthreads; i++) {
       threadData[i].gval = gval; 
       threadData[i].nthreads = nthreads;
       threadData[i].img = img;
       threadData[i].outRef = outRef;
    }

    
    printf("Command: %s\n", argv[0]);
    printf("Filtering image '%s' using attribute area with lambda=%f.\n", imgfname, lambda);
    printf("FreeImage version: %s\n", FreeImage_GetVersion());
    printf("Image: Width=%ld Height=%ld Depth=%ld Size=%ld Size2D=%ld. ", img.width, img.height, img.depth, img.size, img.size2D);
    printf("nthreads for sorting and quantize the max tree: %d.\n", nthreads);
    printf("Size of int=%lu. Size of long=%ld. Size of float=%ld. Size of double=%ld.\n", sizeof(int), sizeof(long), sizeof(float), sizeof(double));
    
    //hmin = FindMin(gval); hmax = FindMax(gval);
    //printf("Min=%d. Max=%d.\n", hmin, hmax);
    
	/**************************************************************/ 
    printf("/*** Sort the pixels ***/\n");  
     
    timings.start_sort = times(&timings.tstruct);    
    // run radix sort
    pixel_t *SORTED = RunRadixSortUnsigned(abs((int)ceil((double)bitsPerPixel/(16))), img, gval, nthreads);	
    for(i = 0; i < nthreads; i++) {
       threadData[i].SORTED = SORTED;
    }
	timings.end_sort = times(&timings.tstruct);  
    
    /******************************************************************/     
    printf("/*** Calculate the quantized image ***/\n");
    timings.start_quimg = times(&timings.tstruct);
    CalculateQuantizedImage(inputQTZLEVELS, threadData);
    timings.end_quimg = times(&timings.tstruct);
    timings.musec = (float)(timings.end_quimg - timings.start_quimg)/((float)timings.tickspersec);
		
	/******************************************************************/ 	
        
    printf("/*** Build the max tree of the quantized image. (threads %d)***/\n", nthreads);    
    timings.start_qutree = times(&timings.tstruct);
    BuildMaxTreeOfQuantizedImage(threadData);
    timings.end_qutree = times(&timings.tstruct);   
    int nthreadsRef = threadData->numQTZLEVELS;   
    printf("Pilot max-tree built. (numQTZLEVELS=%d).\n", threadData->numQTZLEVELS);
    
    /******************************************************************/ 	
    printf("/*** Refinement phase (threads %d) ***/\n", nthreadsRef);

        timings.start_ref = times(&timings.tstruct);	
        RefineTreeBerger(nthreadsRef, threadData);
        timings.end_ref = times(&timings.tstruct);

	printf("Refined max-tree built.\n"); 
	   
    /******************************************************************/
	printf("Init filtering\n");   
    timings.start_filt = times(&timings.tstruct);
    Filter(threadData, lambda);      
    timings.end_filt = times(&timings.tstruct);
    printf("End filtering\n");
    /******************************************************************/ 	
    printTimings(timings);
	
	/**** Write output image   ****/   

    WriteFITS3D(outfname, imgfname, threadData->outRef, img);
    printf("Image written to '%s'\n", outfname);
	
    printf("\n");
    FreeThreadData(threadData);

    FreeImage_DeInitialise();
    return(0);
}
