/* It builds a max tree of the quantized image + support functions for
 * threads: Barrier, semaphores, ...*/

#ifndef _MAXTREE_H
#define _MAXTREE_H

#define QueueFirst(hq,h)  (hq[h].Pixels[hq[h].Head++])
#define QueueAdd(hq,h,p)  hq[h].Pixels[hq[h].Tail++] = p
#define QueueNotEmpty(hq,h)  (hq[h].Head != hq[h].Tail)



void *SafeMalloc(int n);
void *SafeCalloc(int nmemb, pixel_t size);
void Psa(int p);
void Vsa(int p);
void Barrier(int self, int nthreads);
pixel_t levroot(pixel_t x, int *gvalues);
pixel_t Par(pixel_t x, int *gvalues);
void levrootfix(pixel_t lwb, pixel_t upb, int *gvalues);
 
Queue *QueueCreate(pixel_t imgsize, int numQTZLEVELS);
void FillPixelsInQueue(Queue *hq, pixel_t *numpixelsperlevel, int numQTZLEVELS);
void QueueDelete(Queue *hq);

int GetNeighbors(pixel_t p, pixel_t x, pixel_t y, pixel_t z, pixel_t *neighbors, pixel_t lwb, pixel_t upb, ImageProperties img);
int LocalTreeFlood(int self,  Queue *set, pixel_t *lero, int lev, long *thisarea, MaxNode *nodes, int *gvalues, GeneralThreadData *threadData);
void Connect(pixel_t x, pixel_t y, MaxNode *node, int *gvalues, greyval_t *gval);
void Fuse(int self, int i, MaxNode *node, int *gvalues, ImageProperties img, greyval_t *gval);
void MakeThreadData(int numthreads, GeneralThreadData *data);
void FreeThreadQueues(GeneralThreadData *data, int nthreads);
void BuildQuantizedTree(GeneralThreadData *thdata, int nthreads);
void *ccaf(void *arg);
 
int BuildMaxTreeOfQuantizedImage();

#endif
