#ifndef _QUANTIZEDIMAGE_H
#define _QUANTIZEDIMAGE_H

typedef struct
{
    int self;
} ThreadQntImgData;

int CreateQuantizedImage(int inputQTZLEVELS, GeneralThreadData *threadData);
void CalculateQuantizedImage(int inputQTZLEVELS, GeneralThreadData *threadData);
int WorkOutBoundaries(int inputQTZLEVELS, GeneralThreadData *threadData);

void CreateQuantizedImageInParallel(int numQtzLevsUsed, GeneralThreadData *threadData);
void MakeThreadPixelPositions(int numthreads, GeneralThreadData *data);
void MakeThreadQntImgData(int numthreads, GeneralThreadData *threadData);
void RunThreadsWritingQuantizedImage(GeneralThreadData *thdata, int nthreads);
void *wqi(void *arg);

#endif
