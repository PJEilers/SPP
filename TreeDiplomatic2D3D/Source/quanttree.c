#include "TreeDiplomatic.h"

pthread_mutex_t samut[MAXTHREADS];
pthread_cond_t  sacv[MAXTHREADS];
int             saval[MAXTHREADS];
pthread_mutex_t barriermutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t barriercv = PTHREAD_COND_INITIALIZER;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
int barcnt = 0;

/************************** safe malloc and calloc ************************************/



void *SafeMalloc(int n)
{
    void *ptr;
    pthread_mutex_lock(&mutex);
    ptr = malloc(n);
    if (ptr==NULL)
    {
        fprintf (stderr, "Error: out of memory. "
                 "Could not allocate %d bytes\n", n);
    }
    pthread_mutex_unlock(&mutex);
    return ptr;
}

void *SafeCalloc(int nmemb, pixel_t size)
{
    void *ptr;
    pthread_mutex_lock(&mutex);
    ptr = calloc(nmemb, size);
    if (ptr==NULL)
    {
        fprintf (stderr, "Error: out of memory. Could not allocate %ld bytes\n", nmemb*size);
    }
    pthread_mutex_unlock(&mutex);
    return ptr;
}

/************************** semaphore ************************************/

void Psa(int p)
{
    pthread_mutex_lock(&samut[p]);
    while (saval[p] <= 0)
        pthread_cond_wait(&sacv[p], &samut[p]);
    saval[p]--;
    pthread_mutex_unlock(&samut[p]);
}

void Vsa(int p)
{
    pthread_mutex_lock(&samut[p]);
    saval[p]++;
    pthread_mutex_unlock(&samut[p]);
    pthread_cond_broadcast(&sacv[p]);
}

/************************** barrier ************************************/


void Barrier(int self, int nthreads)
{
    pthread_mutex_lock(&barriermutex);
    barcnt++;
    if (barcnt == nthreads)
    {
        barcnt = 0;  /* for reuse of routine */
        pthread_cond_broadcast (&barriercv);
    }
    else
    {
        pthread_cond_wait (&barriercv, &barriermutex);
    }
    pthread_mutex_unlock(&barriermutex);
}

/************************** level root ************************************/

pixel_t levroot(pixel_t x, int *gvalues, MaxNode *node_qu)
{
    pixel_t r=x, y ;
    greyval_t gv=gvalues[x];
    if (r==BOTTOM) return BOTTOM;
    while ((node_qu[r].parent!=BOTTOM) && (gv==gvalues[node_qu[r].parent]))
        r = node_qu[r].parent;

    while (x!=r)
    {
        y=node_qu[x].parent;
        node_qu[x].parent=r;
        x=y;
    }
    return r;
}

pixel_t Par(pixel_t x, int *gvalues, MaxNode *node_qu)
{
    return (pixel_t) levroot(node_qu[x].parent, gvalues, node_qu);
}

void levrootfix(pixel_t lwb, pixel_t upb, int *gvalues, MaxNode *node_qu)
{
    pixel_t z, u, x;
    
    for (x=lwb; x<upb; x++)
    {
        u = levroot(x, gvalues, node_qu);
        if (x!=u) node_qu[x].parent=u;
        else
        {
            z = Par(x, gvalues, node_qu);
            node_qu[x].parent=z;
        }
    }
}

/************************** queue ************************************/
Queue *QueueCreate(pixel_t imgsize, int numQTZLEVELS)
{
    Queue *hq;

    hq = calloc(numQTZLEVELS, sizeof(Queue));
    if (hq==NULL)
    {
        fprintf (stderr, "Out of memory!");
        return(NULL);
    }
    hq->Pixels = calloc(imgsize, sizeof(pixel_t));
    if (hq->Pixels==NULL)
    {
        free(hq);
        fprintf (stderr, "Out of memory!");
        return(NULL);
    }
    return(hq);
} /* QueueCreate */

void FillPixelsInQueue(Queue *hq, pixel_t *numpixelsperlevel, int numQTZLEVELS)
{
    int i;
    hq->Head = hq->Tail = 0;

    for (i=1; i<numQTZLEVELS; i++)
    {
        hq[i].Pixels = hq[i-1].Pixels + numpixelsperlevel[i-1];
        hq[i].Head = hq[i].Tail = 0;
    }
} /* FillPixelsInQueue */

void QueueDelete(Queue *hq)
{
    free(hq->Pixels);
    free(hq);
} /* QueueDelete */

/****************** Construction of local tree using algorithm flood****************************/
int GetNeighbors(pixel_t p, pixel_t x, pixel_t y, pixel_t z, 
		 pixel_t *neighbors, pixel_t lwb, pixel_t upb, ImageProperties img)
{
   int n=0;
   long width = img.width;
   long size2D = img.size2D;
   long height = img.height;
   long depth = img.depth;
   x = p % width;
   if (x<(width-1))       neighbors[n++] = p+1;
   if (y>0)               neighbors[n++] = p-width; // never exec in 2D
   if (x>0)               neighbors[n++] = p-1;
   if (y<(height-1))      neighbors[n++] = p+width; // never exec in 2D
   if (depth>1) 
   {
      if (z>lwb)          neighbors[n++] = p-size2D;
      if (z<(upb-1))           neighbors[n++] = p+size2D;
   }
   return(n);
}

/* parameter MaxNode *nodes must be the par array of the quantized image: 'node_qu';
 * parameter greyval_t *gvalues must be the quantized image: 'gval_qu' */
int LocalTreeFlood(int self,  Queue *set, pixel_t *lero, int lev, void **thisattr, MaxNode *nodes, int *gvalues, GeneralThreadData *threadData)
{

    pixel_t neighbors[CONNECTIVITY];
    int numQTZLEVELS = threadData->numQTZLEVELS;
    int nthreads = threadData->nthreads;
    long width = threadData->img.width;
    long size2D = threadData->img.size2D;
    bool *reached_qu = threadData->reached_qu;
    greyval_t *gval = threadData->gval;
    AuxDataStore *store = &threadData->store;
    
    void *attr = NULL, *childattr;
    
    pixel_t p, q, x, y, z, lwb, upb;
    int fq;
    int numneighbors, i, m;

    lwb = LWB(self, nthreads, threadData->img.size2D, threadData->img.depth);
    upb = UPB(self, nthreads, threadData->img.size2D, threadData->img.depth);
    while(QueueNotEmpty(set, lev))
    {
        p = QueueFirst(set, lev);
        
        x = p % width; 
        y = (p % size2D) / width;
        z = p / size2D;
        numneighbors = GetNeighbors(p, x, y, z, neighbors, lwb/size2D, upb/size2D, threadData->img);
        
        if (attr)  {
            AddToAuxData(attr, x, y, z);
        } else {
            attr = NewAuxData(store, x, y, z);
            if (attr==NULL)  return(NUMBUCKETS);
            if (*thisattr) MergeAuxData(attr, *thisattr);
        }
        
        for (i=0; i<numneighbors; i++)
        {
            q = neighbors[i];
            if (!reached_qu[q])
            {
                fq = gvalues[q];
                reached_qu[q] = true;

                if (lero[fq]==BOTTOM)
		  {
		    lero[fq] = q;
		  }
                else 
		  {
		    if( ( (gval[q] < gval[lero[fq]]) || ( (gval[q] == gval[lero[fq]]) && (q < lero[fq]) )) )
		      {
			nodes[lero[fq]].parent = q; //new lero
			lero[fq] = q; //new lero						
		      }
		    nodes[q].parent = lero[fq];
		  }
					
                QueueAdd(set, fq, q);
                
                if (fq > lev)
                {
                    childattr = NULL;
                    do
                    {
                        fq = LocalTreeFlood(self, set, lero, fq, &childattr, nodes, gvalues, threadData);
                        if(fq>=numQTZLEVELS)
                        {
                            DeleteAuxData(attr);
                            return(fq);
                        }
                    }
                    while (fq!=lev);
                    MergeAuxData(attr, childattr);
                }
            }
        }
    }        
    m = lev-1;
    while ((m>=0) && (lero[m]==BOTTOM))
    {
		m--;
	}

    if (m>=0) nodes[lero[lev]].parent = lero[m];
    // else nodes[lero[lev]].parent = BOTTOM;
    nodes[lero[lev]].Attribute = attr;
    lero[lev] = BOTTOM;
    *thisattr = attr;
    return(m);
} /* MaxTreeFlood */

/************************ Merge trees *************************************/

void Connect(pixel_t x, pixel_t y, MaxNode *node, int *gvalues, greyval_t *gval, AuxDataStore *store)
{

    
    void *cor = NULL, *copa = NULL;
    
    pixel_t h, z;

    x = levroot(x, gvalues, node);
    y = levroot(y, gvalues, node);
    //if (gvalues[x] < gvalues[y]) // changed to keep the levelroot of each component corresponding to the pixel with the minimum intensity value in the original image.
    if ( gval[x] < gval[y] || ((gval[x] == gval[y]) && (x < y)) )
    {
        h=x;
        x=y;
        y=h;
    }
    while ((x!=y) && (y!=BOTTOM))
    {
        z = Par(x, gvalues, node);
        //if ((z!=BOTTOM) && (gvalues[z]>=gvalues[y]))
        //if ( (z!=BOTTOM) && ((gval[z]>=gval[y])))
        if ( (z!=BOTTOM) && ( (gval[z]>gval[y]) || ((gval[z]==gval[y]) && (z>=y)) ) ) //z>=y
        {
            if(cor) MergeAuxData(node[x].Attribute, cor);
            x = z;
        }
        else
        {
            
            if(cor) MergeToAuxData(store, &copa, node[x].Attribute, cor);
            else CloneAuxData(store, &copa, node[x].Attribute);
            CloneAuxData(store, &cor, node[x].Attribute);
            CloneAuxData(store, &node[x].Attribute, copa);
                   
            
            node[x].parent = y;
            x = y;
            y = z;
        }
    }
    if (y==BOTTOM)
    {
        while(x!=BOTTOM)
        {
            if (cor) MergeAuxData(node[x].Attribute, cor); 
            x = Par(x, gvalues, node);
        }
    }
    
    if (cor) DeleteAuxData(cor);
    if (copa) DeleteAuxData(copa);
}


void Fuse(int self, int i, MaxNode *node, int *gvalues, ImageProperties img, greyval_t *gval, int nthreads, AuxDataStore *store) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
    pixel_t p, q, x, y;

    greyval_t prevmin, curmin;
    long width = img.width;
    long height = img.height;
    long size2D = img.size2D;

    /* get the horizontal boundary */

    p  = LWB(self+i, nthreads, size2D, img.depth);
	q = p - size2D;

    x = p % width;
    y = (p /width) % height;

    /*  printf("Region %d merger with %d: (%d,%d,%d)\n",self, self+i, x,y,z);*/

    for ( y = 0 ; y < height ; y++ )
    {
        Connect(p, q, node, gvalues, gval, store);
        //prevmin = MIN(gvalues[p], gvalues[q]);
        prevmin = MIN(gval[p], gval[q]);
        
        p++;
        q++;
        for ( x = 1 ; x < width ; x++, p++, q++ )
        {
            //curmin = MIN(gvalues[p], gvalues[q]); // changed from GVALUES to GVAL to keep the levelroot of each component corresponding to the pixel with the minimum intensity value in the original image.
            curmin = MIN(gval[p], gval[q]);
			if (curmin > prevmin)
            {
                Connect(p, q, node, gvalues, gval, store);
            }
            prevmin = curmin;
        }
    }
}

void MakeThreadData(int numthreads, GeneralThreadData *data)
{
    int i;

    for (i=0; i<numthreads; i++)
    {
        data[i].thisqueue= QueueCreate( UPB(i, numthreads, data->img.size2D, data->img.depth)-LWB(i, numthreads, data->img.size2D, data->img.depth), data->numQTZLEVELS);
    }
}

void FreeThreadQueues(GeneralThreadData *data, int numthreads)
{
    int i;
    for (i=0; i<numthreads; i++)
        QueueDelete(data[i].thisqueue);
}


void *ccaf(void *arg)
{
    GeneralThreadData *thdata = (GeneralThreadData *) arg;
    int numQTZLEVELS = thdata-> numQTZLEVELS;
    bool *reached_qu = thdata->reached_qu;
    MaxNode *node_qu = thdata->node_qu;
    int self = thdata->self, q, i;
    int *gval_qu = thdata->gval_qu;
    int nthreads = thdata->nthreads;
    greyval_t *gval = thdata->gval;
    pixel_t x;
    void *attr = NULL;
    pixel_t numpixelsperlevel[numQTZLEVELS], lero[numQTZLEVELS], xm;
    Queue *set = thdata->thisqueue;
    AuxDataStore *store = &thdata->store;
    
    for (i=0; i<numQTZLEVELS; i++)
    {
        numpixelsperlevel[i] = 0;
        lero[i] = BOTTOM;
    }

    xm = LWB(self, nthreads, thdata->img.size2D, thdata->img.depth);

    for (x=xm; x<UPB(self, nthreads, thdata->img.size2D, thdata->img.depth); x++)
    {
        numpixelsperlevel[gval_qu[x]]++;
        node_qu[x].parent = BOTTOM;
        //if (gval_qu[xm]>gval_qu[x]) xm = x;
        if (gval[xm]>gval[x]) xm = x; 
    }
    FillPixelsInQueue(set, numpixelsperlevel, numQTZLEVELS);

    QueueAdd(set, gval_qu[xm], xm);
    reached_qu[xm] = true;
    
    // the level root for gval_qu[] for this thread is equals to 
    // the minimum original value for the quantized component.
    lero[gval_qu[xm]] = xm;

    LocalTreeFlood(self, set, lero, gval_qu[xm], &attr, node_qu, gval_qu, thdata);
    /*for (i=0; i<size; i++)
      if(node_qu[i].parent == BOTTOM)
	printf("%d %ld %ld\n", i ,node_qu[i].parent, SORTED[i]);*/
    i = 1;
    q = self;
    while ((self+i<nthreads) && (q%2 == 0))
    {
        Psa(self+i);  /* wait to glue with righthand neighbor */
        Fuse(self, i, node_qu, gval_qu, thdata->img, gval, nthreads, store);
        i = 2*i;
        q = q/2;
    }
    if (self != 0)
    {
        Vsa(self);  /* signal lefthand neighbor */
    }
    return NULL;
}

void BuildQuantizedTree(GeneralThreadData *thdata, int nthreads)
{
    pthread_t threadID[MAXTHREADS];
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, ccaf, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

int BuildMaxTreeOfQuantizedImage(GeneralThreadData *threadData)
{
    // Build the max tree of the quantized image
    int attrib = threadData->attrib;
    AttribStruct *Attribs = getAttribs();
    
    NewAuxData = Attribs[attrib].NewData;
    DeleteAuxData = Attribs[attrib].DeleteData;
    AddToAuxData = Attribs[attrib].AddToData;  
    MergeAuxData = Attribs[attrib].MergeData;
    MergeToAuxData = Attribs[attrib].MergeToData;
    CloneAuxData = Attribs[attrib].CloneData;
    
    
    int nthreads = threadData->nthreads;
    int i;
    long size = threadData->img.size;
    int *gval_qu = threadData->gval_qu;
    greyval_t *gval = threadData->gval;
    MaxNode *node_qu = calloc((size_t)size, sizeof(MaxNode));
    bool *reached_qu = calloc((size_t)size, sizeof(bool));
    if (node_qu==NULL)
    {
        fprintf(stderr, "out of memory! \n");
        free(gval);
        free(gval_qu);
        return(-1);
    }
    if (reached_qu==NULL)
    {
        fprintf(stderr, "out of memory!\n");
        free(node_qu);
        free(gval_qu);
        return(-1);
    }

    
    for(i = 0; i < nthreads; i++) {
        threadData[i].reached_qu = reached_qu;
        threadData[i].node_qu = node_qu;
    }
    

    MakeThreadData(nthreads, threadData);
    BuildQuantizedTree(threadData, nthreads);
    FreeThreadQueues(threadData, nthreads);  
    return 0;
}
