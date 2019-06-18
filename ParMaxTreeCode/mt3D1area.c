#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

#include "ReadWriteData.h"
#include "Volume.h"
#include "avs_io.h"

#define PI 3.14159265358979323846
#define false 0
#define true  1
#define MAXTHREADS 128
#define MAXLEVELS  4096          /* altered to fit in 12 bit data */
#define CONNECTIVITY  6

/* typedef short bool;
   typedef unsigned int uint; 
   typedef unsigned long ulong; */
typedef unsigned char ubyte;

typedef ushort voxel;

#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))
#define LWB(self) (size2D*(((self)*depth)/nthreads))
#define UPB(self) (size2D*(((self+1)*depth)/nthreads))
#define bottom (-1)

int nthreads;
pthread_t threadID[MAXTHREADS];

ubyte imgftype;
ulong width, height, depth, size;  /* precondition: width <= size/nthreads */
ulong size2D;
voxel *gval=NULL, *out=NULL;
VolumeStruct VolStruct;
int attrib, decision=1; /*  the original default value of decision is 3; */
double lambda;

typedef struct MaxNode 
{ 
   ulong parent;
   ulong Area;
   int filter; /* gray level after filtering */
} MaxNode;

MaxNode *node; 
bool *reached;
bool *remaining; 

/************************** safe malloc and calloc ************************************/

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void *SafeMalloc(int n) {
  void *ptr;
  pthread_mutex_lock(&mutex);
  ptr = malloc(n);
  if (ptr==NULL) {
    fprintf (stderr, "Error: out of memory. "
	     "Could not allocate %d bytes\n", n);
  }
  pthread_mutex_unlock(&mutex);
  return ptr;
}

void *SafeCalloc(int nmemb, int size) {
  void *ptr;
  pthread_mutex_lock(&mutex);
  ptr = calloc(nmemb, size);
  if (ptr==NULL) {
    fprintf (stderr, "Error: out of memory. Could not "
	     "allocate %d bytes\n", nmemb*size);
  }
  pthread_mutex_unlock(&mutex);
  return ptr;
}

/************************** semaphore ************************************/

pthread_mutex_t samut[MAXTHREADS];
pthread_cond_t  sacv[MAXTHREADS];
int             saval[MAXTHREADS];

void Psa(int p) {
  pthread_mutex_lock(&samut[p]);
  while (saval[p] <= 0)
    pthread_cond_wait(&sacv[p], &samut[p]);
  saval[p]--;
  pthread_mutex_unlock(&samut[p]);
}

void Vsa(int p) {
  pthread_mutex_lock(&samut[p]);
  saval[p]++;
  pthread_mutex_unlock(&samut[p]);
  pthread_cond_broadcast(&sacv[p]);
}

/************************** barrier ************************************/

pthread_mutex_t barriermutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t barriercv = PTHREAD_COND_INITIALIZER;
int barcnt = 0;

void Barrier(int self) {
  pthread_mutex_lock(&barriermutex);
  barcnt++;
  if (barcnt == nthreads) {
    barcnt = 0;  /* for reuse of routine */
    pthread_cond_broadcast (&barriercv);  
  } else {
    pthread_cond_wait (&barriercv, &barriermutex);
  }
  pthread_mutex_unlock(&barriermutex);
}

/************************** level root ************************************/

ulong levroot(ulong x) {
  int r=x, y ,gv=gval[x];
  if (r==bottom) return bottom;
  while ((node[r].parent!=bottom) && (gv==gval[node[r].parent]))
      r = node[r].parent;

  while (x!=r){
    y=node[x].parent;
    node[x].parent=r;
    x=y;
  }
  return r;
}

ulong Par(ulong x) {
  return levroot(node[x].parent);
}

void levrootfix(ulong lwb, ulong upb) {
  ulong z, u, x;

  for (x=lwb; x<upb;x++) {
    u = levroot(x);
    if (x!=u) node[x].parent=u;
    else {
      z = Par(x); 
      node[x].parent=z;
    }
  }
}

void MaxTreeAreaFilter(int self, double lambda)
{
   ulong v, u, w, parent, lwb, upb;
   int val;
 
   lwb = LWB(self); upb = UPB(self);
   for (v=lwb; v<upb; v++) {
      if (node[v].filter == -1) {
         w = v;
         parent = node[w].parent;
         while ((parent != bottom) && (node[w].filter == -1) && 
            ((gval[w] == gval[parent]) || (node[w].Area < lambda))) {
                  w = parent;
                  parent = node[w].parent;
         }
         if (node[w].filter != -1) val = node[w].filter;
         else if (node[w].Area >= lambda) val = gval[w]; 
         else val = 0; /* criterion cannot be satisfied */
         u = v;
         while (u!=w) {
           if ((lwb<=u) && (u<upb)) node[u].filter = val;
           u = node[u].parent;
           }
         if ((lwb<=w) && (w<upb)) node[w].filter = val;
       }
       out[v] = node[v].filter;
   }
} /* MaxTreeAreaFilter */



/************************** queue ************************************/

typedef struct Queue
{
   ulong *Pixels;
   ulong Head;
   ulong Tail; /* First free place in queue, empty if Head=Tail */
} Queue;

Queue *QueueCreate(ulong imgsize)
{
   Queue *hq;

   hq = calloc(MAXLEVELS, sizeof(Queue));
   if (hq==NULL) {
       fprintf (stderr, "Out of memory!");
       return(NULL);
   }
   hq->Pixels = calloc(imgsize, sizeof(ulong));
   if (hq->Pixels==NULL)
   {
      free(hq);
      fprintf (stderr, "Out of memory!");
      return(NULL);
   }
  return(hq);
} /* QueueCreate */

void FillPixelsInQueue(Queue *hq, ulong *numpixelsperlevel)
{
  int i;
  hq->Head = hq->Tail = 0;

  for (i=1; i<MAXLEVELS; i++)
    {
      hq[i].Pixels = hq[i-1].Pixels + numpixelsperlevel[i-1];
      hq[i].Head = hq[i].Tail = 0;
    }
} /* FillPixelsInQueue */

#define QueueFirst(hq,h)  (hq[h].Pixels[hq[h].Head++])
#define QueueAdd(hq,h,p)  hq[h].Pixels[hq[h].Tail++] = p
#define QueueNotEmpty(hq,h)  (hq[h].Head != hq[h].Tail)

void QueueDelete(Queue *hq)
{
   free(hq->Pixels);
   free(hq);
} /* QueueDelete */

/****** Image read/write functions ******************************/

short ImagePGMAsciiRead(char *fname)
{
   FILE *infile;
   ulong i;
   int c;

   infile = fopen(fname, "r");
   if (infile==NULL) {
      fprintf (stderr, "Error: Can't read the ASCII file: %s !", fname);
      return(0);
   }
   fscanf(infile, "P2\n");
   while ((c=fgetc(infile)) == '#')
      while ((c=fgetc(infile)) != '\n');
   ungetc(c, infile);
   fscanf(infile, "%lu %lu\n255\n", &width, &height);
   size = width*height;
   gval = malloc(size*sizeof(voxel));
   if (gval==NULL) {
      fprintf (stderr, "Out of memory!");
      fclose(infile);
      return(0);
   }
   for (i=0; i<size; i++)
   {
      fscanf(infile, "%d", &c);
      gval[i] = c;
   }
   fclose(infile);
   return(1);
} /* ImagePGMAsciiRead */


short ImagePGMBinRead(char *fname)
{
   FILE *infile;
   int c, i;
   ubyte *buf = NULL;

   infile = fopen(fname, "rb");
   if (infile==NULL) {
      fprintf (stderr, "Error: Can't read the binary file: %s !", fname);
      return(0);
   }
   fscanf(infile, "P5\n");
   while ((c=fgetc(infile)) == '#')
      while ((c=fgetc(infile)) != '\n');
   ungetc(c, infile);
   fscanf(infile, "%lu %lu\n255\n", &width, &height);
   size = width*height;
   buf = malloc(size);
   if (buf==NULL) {
     fprintf (stderr, "Out of memory!");
     fclose(infile);
     return(0);
   } 
   gval = malloc(size*sizeof(voxel));
   if (gval==NULL) {
     fprintf (stderr, "Out of memory!");
     fclose(infile);
     return(0);
   } 
   fread(buf, 1, size, infile);

   for (i=0; i<size; i++)
     gval[i]=buf[i];
   free(buf);
   fclose(infile);
   return(1);
} /* ImagePGMBinRead */


short ImagePGMRead(char *fname)
{
   FILE *infile;
   char id[4];

   infile = fopen(fname, "r");
   if (infile==NULL) {
      fprintf (stderr, "Error: Can't read the image: %s !", fname);
      return(0);
   }
   fscanf(infile, "%3s", id);
   fclose(infile);
   if (strcmp(id, "P2")==0) return(ImagePGMAsciiRead(fname));
   else if (strcmp(id, "P5")==0) return(ImagePGMBinRead(fname));
   else {
     fprintf (stderr, "Unknown type of the image!");
     return(0);
   }
} /* ImagePGMRead */

int ImagePGMBinWrite(char *fname)
{
   FILE *outfile;
   ubyte *buf=NULL;
   int i;

   outfile = fopen(fname, "wb");
   if (outfile==NULL) {
      fprintf (stderr, "Error: Can't write the image: %s !", fname);
      return(-1);
   }
   fprintf(outfile, "P5\n%ld %ld\n255\n", width, height);

   buf = malloc(size);
   if (buf==NULL) {
     fprintf (stderr, "Out of memory!");
     fclose(outfile);
     return(-1);
   } 
   for (i=0;i<size;i++)
     buf[i]=out[i];
   fwrite(buf, 1, (size_t)(size), outfile);
   free(buf);
   fclose(outfile);
   return(0);
} /* ImagePGMBinWrite */


/****************** Construction of local tree using algorithm flood****************************/

int GetNeighbors(ulong p, int x, int y, int z, 
		 ulong *neighbors, ulong lwb, ulong upb)
{

   int n=0;

   x = p % width;

   if (x<(width-1))       neighbors[n++] = p+1;
   if (y>0)               neighbors[n++] = p-width;
   if (x>0)               neighbors[n++] = p-1;
   if (y<(height-1))      neighbors[n++] = p+width;
   if (depth>1) {
      if (z>lwb)          neighbors[n++] = p-size2D;
      if (z<(upb-1))           neighbors[n++] = p+size2D;
   }
   return(n);
} /* GetNeighbors */

int LocalTreeFlood(int self,  Queue *set, ulong *lero, int lev, ulong *thisarea)
{
   ulong neighbors[CONNECTIVITY];
   ulong area = *thisarea, childarea;
   ulong p, q, x, y, z, lwb, upb;
   int fq;
   int numneighbors, i, m;

   lwb = LWB(self);
   upb = UPB(self);

   while(QueueNotEmpty(set, lev)) {
      area++;
      p = QueueFirst(set, lev);
      x = p % width; 
      y = (p % size2D) / width;
      z = p / size2D;
      numneighbors = GetNeighbors(p, x, y, z, neighbors, lwb/size2D, upb/size2D);
      for (i=0; i<numneighbors; i++) {
         q = neighbors[i];
         if (!reached[q]) {
            fq = gval[q];            
            reached[q] = true;
            if (lero[fq]==bottom) lero[fq] = q;
            else node[q].parent = lero[fq];
            QueueAdd(set, fq, q);
            if (fq > lev) {
               childarea = 0;
               do {
                  fq = LocalTreeFlood(self, set, lero, fq, &childarea);
                  if (fq>=MAXLEVELS) {
                     return(fq);
                  }
               } while (fq!=lev);
               area += childarea;
            }
         }
      }
   }
   m = lev-1;
   while ((m>=0) && (lero[m]==bottom))  m--;
   if (m>=0) node[lero[lev]].parent = lero[m];
   node[lero[lev]].Area = area;
   lero[lev] = bottom;
   *thisarea = area;
   return(m);
} /* MaxTreeFlood */

/************************ Merge trees *************************************/

void Connect(ulong x, ulong y) {

  ulong area = 0, area1 = 0;
  ulong h, z;

  x = levroot(x);
  y = levroot(y);
  if (gval[x] < gval[y]) {
    h=x; x=y; y=h;
  }
  while ((x!=y) && (y!=bottom)) {
    z = Par(x);
    if ((z!=bottom) && (gval[z]>=gval[y])) {
      node[x].Area += area;
      x = z;
    } else {
      area1 = node[x].Area + area;
      area = node[x].Area ;
      node[x].Area =area1;
      node[x].parent = y;
      x = y;
      y = z;
    }
  }
  if (y==bottom) {
    while(x!=bottom) {
      node[x].Area += area;
      x = Par(x);
    }
  }
}

void Fuse2(int self, int i) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
  ulong lwb, upb;
  ulong p, xm, x, y;
   
  lwb = LWB(self);
  upb = UPB(self+2*i-1); 
  if (upb>size) upb=size;

  /* get the horizontal boundary */
  xm = LWB(self+i); 

  x = xm;
  p = x % width;
  if ((p>0) && (x-1>=lwb)) {
     y = x-1;
     Connect(x, y);
  }  

  for (p=0; p<width && x<upb; p++){
      if (x>=lwb+width) {
         y= x-width;
         Connect(x, y);
      }
      x++;
  }

  if (depth>1) {
    x = xm;
    for (p=0; p<size2D && x<upb; p++){
      if (x>=lwb+size2D) {
	y= x-size2D;
	Connect(x, y);
      }
      x++;
    }
  }

  /* levrootfix(lwb, upb); */
}


void Fuse(int self, int i) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
  ulong p, q, x, y, z;

  voxel prevmin, curmin;

  /* get the horizontal boundary */ 

  p  = LWB(self+i);
  q = p - size2D;

  x = p % width;
  y = (p /width) % height;
  z = p /size2D;
  
  /*  printf("Region %d merger with %d: (%d,%d,%d)\n",self, self+i, x,y,z);*/

  for ( y = 0 ; y < height ; y++ ){
    Connect(p, q);
    prevmin = MIN(gval[p],gval[q]);
    p++;
    q++;
    for ( x = 1 ; x < width ; x++, p++, q++ ){
      curmin = MIN(gval[p],gval[q]);
      if (curmin > prevmin){
	Connect(p, q);
      }
      prevmin = curmin;    
    }  
  }

  /* levrootfix(lwb, upb); */
}

/****************** Concurrent construction and filter of Maxtree  ****************************/

typedef struct { 
  int self;
  Queue *thisqueue;
} ThreadData;

ThreadData *MakeThreadData(int numthreads){
  ThreadData *data = malloc(numthreads *sizeof(ThreadData));
  int i;

  for (i=0; i<numthreads; i++){
    data[i].self=i;
    data[i].thisqueue= QueueCreate( UPB(i)-LWB(i));
    printf("Thread %d, LWB = %d, UPB = %d\n",i,LWB(i),UPB(i));
  }    
  return(data);
}

void FreeThreadData(ThreadData *data, int numthreads)
{
  int i;
  for (i=0;i<numthreads; i++)
    QueueDelete(data[i].thisqueue);
  free(data);
}

void *ccaf(void *arg) {
  ThreadData *thdata = (ThreadData *) arg;
  int self = thdata->self, q, i;
  ulong x, area=0;
  ulong numpixelsperlevel[MAXLEVELS], lero[MAXLEVELS], xm;
  Queue *set = thdata->thisqueue;

  for (i=0; i<MAXLEVELS; i++) {
      numpixelsperlevel[i] = 0;
      lero[i] = bottom;
  }

  xm = LWB(self);

  for (x=xm; x<UPB(self); x++) { 
      numpixelsperlevel[gval[x]]++; 
      node[x].parent = bottom;
      node[x].filter = -1;
      if (gval[xm]>gval[x]) xm = x;
  }
  FillPixelsInQueue(set, numpixelsperlevel);

  QueueAdd(set, gval[xm], xm);
  reached[xm] = true;
  lero[gval[xm]] = xm;             

  LocalTreeFlood(self, set, lero, gval[xm], &area);



  i = 1;
  q = self;
  while ((self+i<nthreads) && (q%2 == 0)) {
    Psa(self+i);  /* wait to glue with righthand neighbor */
    Fuse(self, i); 
    i = 2*i;
    q = q/2;
  }
  if (self != 0) {
    Vsa(self);  /* signal lefthand neighbor */
  }
  Barrier(self);
 
  MaxTreeAreaFilter(self, lambda);

  if (decision==2) {   /* thread 0 needs to do some remaining work for decision Max if necessary */
     if (self != 0) {
        Vsa(self);  /* signal thread 0*/
     } else {
       for (i=1; i<nthreads; i++) Psa(i);   /* waiting for all other threads */
       for (x=0; x<size; x++) {    
	 if (remaining[x]) {
            xm = x;
            while (xm!=bottom && node[xm].filter!=gval[xm]) {
	      node[xm].filter = gval[xm];
              out[xm] = gval[xm];
              xm = node[xm].parent;
            }
	 } 
       }           
     }
  }
  return NULL;
}

void BuildTreeAndFilter(ThreadData *thdata, int nthreads) {
  int thread;
  for (thread=0; thread<nthreads; ++thread) {
    pthread_create(threadID+thread, NULL, ccaf, (void *) (thdata + thread));
  }
  for (thread=0; thread<nthreads; ++thread) {
    pthread_join(threadID[thread], NULL);
  }
}

int main (int argc, char *argv[]) {
  
   char *imgfname, *outfname = "out.pgm";
   avs_header avs_head;
   ThreadData *thdata;
   int r;
   ulong i;
   clock_t start;
   struct tms tstruct;
   long tickspersec = sysconf(_SC_CLK_TCK);  
   float musec;

   if (argc<4)
   {
      printf("Usage: %s <nthreads> <input image> <lambda>  [output image] \n", argv[0]);
      exit(0);
   }

   nthreads = MIN(atoi(argv[1]), MAXTHREADS);
   imgfname = argv[2];
   r=0;
   while (imgfname[r]!='.' && imgfname[r]!=0) r++;
   if (imgfname[r+1]=='f' && imgfname[r+2]=='l') imgftype = 3;
   else imgftype = 2;
   lambda = atof(argv[3]);
   if (argc>4)  outfname = argv[4];
   else if (imgftype==3) outfname="out.fld";
   
   if (imgftype==2) {
      if (!ImagePGMRead(imgfname)) return(-1);
      size2D = size;
      depth = 1;
   } else if (readAVS(imgfname, &avs_head, &VolStruct)) {
      if (avs_head.datatype>=4) {
         fprintf (stderr, "The DataType is not byte or short!\n" );
         return(-1);
      }
      width = avs_head.dim1;
      height = avs_head.dim2;
      depth = avs_head.dim3;
      size2D = width*height;
      size = size2D*depth;
      gval = malloc(size*sizeof(voxel));
      if (gval==NULL) {
         fprintf (stderr, "Out of memory!");
         return(-1);
      }
      if (VolStruct.DataType==1)
        for (i=0; i<size; i++) {
	  gval[i] = VolStruct.data.uchar_data[i];
	}
      else
        for (i=0; i<size; i++) {
	  gval[i] = VolStruct.data.ushort_data[i];
	}
	
   }
   printf("Filtering image '%s' using attribute area with lambda=%f\n", imgfname, lambda);
   printf("Image: Width=%ld Height=%ld Depth=%ld\n", width, height, depth);
   printf ("nthreads: %d\n", nthreads);
   
   node = calloc((size_t)size, sizeof(MaxNode));
   if (node==NULL) {
      fprintf(stderr, "out of memory! \n");
      free(gval);
      return(-1);
   }
   reached = calloc((size_t)size, sizeof(bool));
   if (reached==NULL) {
      fprintf(stderr, "out of memory!\n");
      free(node);
      free(gval);
      return(-1);
   }
 
   remaining = calloc((size_t)size, sizeof(bool));
   if (remaining==NULL) {
      fprintf(stderr, "out of memory!\n");
      free(reached);
      free(node);
      free(gval);
      return(-1);
   }
 
   out =  malloc(size*sizeof(voxel));
   if (out==NULL) {
     fprintf(stderr, "Can't create output image! \n");
     free(remaining);
     free(reached);
     free(node);
     free(gval);
     return(-1);
   }
 
   thdata = MakeThreadData(nthreads); 
   printf("Data read, start filtering.\n");
   start = times(&tstruct);
   
   BuildTreeAndFilter(thdata,nthreads);
   
   musec = (float)(times(&tstruct) - start)/((float)tickspersec);

   printf("wall-clock time: %f s\n",musec);
   
   if (imgftype==2) {
     r = ImagePGMBinWrite(outfname);
     free(out);  
   } else {
      if (VolStruct.DataType==1)
        for (i=0; i<size; i++) {
	  VolStruct.data.uchar_data[i]=out[i];
	}
      else
        for (i=0; i<size; i++) {
	  VolStruct.data.ushort_data[i]=out[i];
	}

     r = writeAVS(outfname, &avs_head, &VolStruct);
     DestroyVolume(&VolStruct);
   }
   if (r)  printf("Filtered image written to '%s'\n", outfname);
   
   FreeThreadData(thdata,nthreads);
   free(remaining);
   free(reached);
   free(node);
   free(gval);
   return(0);
} /* main */
