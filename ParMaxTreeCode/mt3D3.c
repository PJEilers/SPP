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
#define MAXTHREADS 64
#define MAXLEVELS  4096           /* altered to fit in 12 bit data */
#define CONNECTIVITY  6

/* typedef short bool;
   typedef unsigned int uint; 
   typedef unsigned long ulong; */
typedef unsigned char ubyte;

typedef ushort voxel;

#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))
#define LWB(self) (((self)*size)/nthreads)
#define UPB(self) (((self+1)*size)/nthreads)
#define bottom (-1)

int nthreads;
pthread_t threadID[MAXTHREADS];

ubyte imgftype;
ulong width, height, depth, size;  /* precondition: width <= size/nthreads */
ulong size2D;
voxel *gval=NULL,  *out=NULL;
VolumeStruct VolStruct;
int attrib, decision=1; /*  the original default value of decision is 3; */
double lambda;

typedef struct MaxNode 
{ 
   ulong parent;
   ulong Area;
   void *Attribute;
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
  if (x==bottom) return bottom;
  while ((node[x].parent!=bottom) && (gval[x]==gval[node[x].parent]))
      x = node[x].parent;
  return x;
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

#include "decision3D.c"

/************************** queue ************************************/

typedef struct Queue
{
   ulong *Pixels;
   int *Neighbors;
   ulong Head;
   ulong Tail; /* First free place in queue, empty if Head=Tail */
} Queue;

Queue *QueueCreate(ulong imgsize, ulong *numpixelsperlevel)
{
   Queue *hq;
   int i;

   hq = SafeCalloc(MAXLEVELS, sizeof(Queue));
   if (hq==NULL) {
       fprintf (stderr, "Out of memory!");
       return(NULL);
   }
   hq->Pixels = SafeCalloc(imgsize, sizeof(ulong));
   if (hq->Pixels==NULL)
   {
      free(hq);
      fprintf (stderr, "Out of memory!");
      return(NULL);
   }
   hq->Neighbors = SafeCalloc(imgsize, sizeof(int));
   if (hq->Neighbors==NULL)
   {
      free(hq->Pixels);
      free(hq);
      fprintf (stderr, "Out of memory!");
      return(NULL);
   }
   hq->Head = hq->Tail = 0;
   for (i=1; i<MAXLEVELS; i++)
   {
      hq[i].Pixels = hq[i-1].Pixels + numpixelsperlevel[i-1];
      hq[i].Neighbors = hq[i-1].Neighbors + numpixelsperlevel[i-1];
      hq[i].Head = hq[i].Tail = 0;
   }
   return(hq);
} /* QueueCreate */

void QueueAdd(Queue *hq, int h, ulong p, int n){
   hq[h].Pixels[hq[h].Tail] = p;
   hq[h].Neighbors[hq[h].Tail] = n;
   hq[h].Tail++;
}

#define QueueNotEmpty(hq,h)  (hq[h].Head != hq[h].Tail)

void QueueDelete(Queue *hq)
{
   free(hq->Pixels);
   free(hq->Neighbors);
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
      return(-2);
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
      return(-1);
   }
   for (i=0; i<size; i++)
   {
      fscanf(infile, "%d", &c);
      gval[i] = c;
   }
   fclose(infile);
   return(0);
} /* ImagePGMAsciiRead */


short ImagePGMBinRead(char *fname)
{
   FILE *infile;
   int c,i;
   ubyte *buf = NULL;

   infile = fopen(fname, "rb");
   if (infile==NULL) {
      fprintf (stderr, "Error: Can't read the binary file: %s !", fname);
      return(-2);
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
     return(-1);
   } 
   fread(buf, 1, size, infile);
   gval = malloc(size * sizeof(voxel));
   if (gval==NULL) {
     fprintf (stderr, "Out of memory!");
     fclose(infile);
     return(-1);
   } 

   for (i=0; i<size; i++)
     gval[i]=buf[i];
   free(buf);
   fclose(infile);
   return(0);
} /* ImagePGMBinRead */


short ImagePGMRead(char *fname)
{
   FILE *infile;
   char id[4];

   infile = fopen(fname, "r");
   if (infile==NULL) {
      fprintf (stderr, "Error: Can't read the image: %s !", fname);
      return(-2);
   }
   fscanf(infile, "%3s", id);
   fclose(infile);
   if (strcmp(id, "P2")==0) return(ImagePGMAsciiRead(fname));
   else if (strcmp(id, "P5")==0) return(ImagePGMBinRead(fname));
   else {
     fprintf (stderr, "Unknown type of the image!");
     return(-3);
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


/****************** Construction of local tree without the recursive call ****************************/

int GetNeighbors(ulong p, ulong *neighbors, ulong lwb, ulong upb)
{
   ulong x;
   int n=0;

   x = p % width;
   if ((x<(width-1)) && (p+1<upb))      neighbors[n++] = p+1;
   if (p>=lwb+width)                    neighbors[n++] = p-width;
   if ((x>0) && (p-1>=lwb))             neighbors[n++] = p-1;
   if (p+width<upb)                     neighbors[n++] = p+width;
   if (depth>1) {
      if (p>=lwb+size2D)          neighbors[n++] = p-size2D;
      if (p+size2D<upb)           neighbors[n++] = p+size2D;
   }
   return(n);
} /* GetNeighbors */

int BuildLocalTree(int self)
{
  ulong neighbors[CONNECTIVITY], numpixelsperlevel[MAXLEVELS], lero[MAXLEVELS];
  ulong p, q, x, y, z, lwb, upb, xm;
  int lev, fq;
  int numneighbors, m, n;
  Queue *wait;

  lwb = LWB(self);
  upb = UPB(self);
  for (m=0; m<MAXLEVELS; m++) {
       numpixelsperlevel[m] = 0;
       lero[m] = bottom;
  }
  xm = lwb;
  for (x=lwb; x<upb; x++) { 
      numpixelsperlevel[gval[x]]++; 
      node[x].parent = bottom;
      node[x].filter = -1;
      if (gval[xm]>gval[x]) xm = x;
  }

  wait = QueueCreate(upb-lwb, numpixelsperlevel);
  if (wait==NULL) return(-1);

  lev = gval[xm];
  QueueAdd(wait, lev, xm, 0);
  reached[xm] = true;
  lero[lev] = xm;             

  while (lev >=0) {
    if (QueueNotEmpty(wait, lev)) { /* Encounter */
        p = wait[lev].Pixels[wait[lev].Head];     /* QueueFirst */
        n = wait[lev].Neighbors[wait[lev].Head];  /* the neighbor that p starts investigating. */
        if (n==0) {
            node[lero[lev]].Area++;
	    x = (p % size2D) % width; 
	    y = (p % size2D) / width;
	    z = p / size2D;
            if (node[lero[lev]].Attribute) AddToAuxData(node[lero[lev]].Attribute, x, y, z); 
            else {
                 node[lero[lev]].Attribute = NewAuxData(x, y, z);
                 if (node[lero[lev]].Attribute==NULL)  return(-1);
     	    }
        } 
        numneighbors = GetNeighbors(p, neighbors, lwb, upb);
        for (m=n; m<numneighbors; m++) {
            q = neighbors[m];
            if (!reached[q]) {
               fq = gval[q];            
               reached[q] = true;
               if (lero[fq]==bottom) lero[fq] = q;
               else node[q].parent = lero[fq];
               QueueAdd(wait, fq, q, 0);
               if (fq > lev) {
 		  wait[lev].Neighbors[wait[lev].Head] = m+1;
                  lev = fq; 
                  break;
               }
            }
        }
        if (m == numneighbors) wait[lev].Head++; /* whether all neighbors of p have been investigated */
     } else { /* Upward */
       m = lev-1;
       while ((m>=0) && (lero[m]==bottom))  m--;
       if (m>=0) {
          node[lero[lev]].parent = lero[m];
          node[lero[m]].Area += node[lero[lev]].Area;
          if (node[lero[m]].Attribute) MergeAuxData(node[lero[m]].Attribute, node[lero[lev]].Attribute); 
          else
             CloneAuxData(&node[lero[m]].Attribute, node[lero[lev]].Attribute); 
       }
       lero[lev] = bottom;
       lev = m;
     }
   } 
  QueueDelete(wait);

  return(0);
} /* BuildLocalTree */ 

/************************ Merge trees *************************************/

void Connect(ulong x, ulong y) {

  ulong area = 0, area1 = 0;
  void *cor = NULL, *copa = NULL;
  ulong h, z;

  x = levroot(x);
  y = levroot(y);
  if (gval[x] < gval[y]) {
    h=x; x=y; y=h;
  }
  while ((x!=y) && (y!=bottom)) {
    z = Par(x);
    if ((z!=bottom) && (gval[z]>=gval[y])) {
      if (cor) MergeAuxData(node[x].Attribute, cor); 
      node[x].Area += area;
      x = z;
    } else {
      if (cor) MergeToAuxData(&copa, node[x].Attribute, cor);
      else CloneAuxData(&copa, node[x].Attribute);
      CloneAuxData(&cor, node[x].Attribute);
      CloneAuxData(&node[x].Attribute, copa);
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
      if (cor) MergeAuxData(node[x].Attribute, cor); 
      node[x].Area += area;
      x = Par(x);
    }
  }
  if (cor) DeleteAuxData(cor);
  if (copa) DeleteAuxData(copa);
}

void Fuse(int self, int i) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
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

/****************** Concurrent construction and filter of Maxtree  ****************************/

void *ccaf(void *arg) {
  int self= (int) arg, q, i;
  ulong xm, x;

  BuildLocalTree(self); 

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
  Decisions[decision].Filter(self, Attribs[attrib].Attribute, lambda);

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

void BuildTreeAndFilter() {
  int thread;
  for (thread=0; thread<nthreads; ++thread) {
    pthread_create(threadID+thread, NULL, ccaf, (void *)thread);
  }
  for (thread=0; thread<nthreads; ++thread) {
    pthread_join(threadID[thread], NULL);
  }
}

int main (int argc, char *argv[]) {
   char *imgfname, *outfname = "out.pgm";
   avs_header avs_head;
   int r;
   ulong i;
   clock_t start;
   struct tms tstruct;
   long tickspersec = sysconf(_SC_CLK_TCK);  
   float musec;

   if (argc<5)
   {
      printf("Usage: %s <nthreads> <input image> <attrib> <lambda> [decision] [output image] \n", argv[0]);
      printf("Where attrib is:\n");
      for (attrib=0; attrib<NUMATTR; attrib++)
      {
         printf("\t%d - %s\n", attrib, Attribs[attrib].Name);
      }
      printf("and decision is:\n");
      for (r=0; r<NUMDECISIONS; r++)
      {
         printf("\t%d - %s", r, Decisions[r].Name);
	 if (r==decision)  printf(" (default)");
	 printf("\n");
      }
      exit(0);
   }

   nthreads = atoi(argv[1]);
   imgfname = argv[2];
   r=0;
   while (imgfname[r]!='.' && imgfname[r]!=0) r++;
   if (imgfname[r+1]=='f' && imgfname[r+2]=='l') imgftype = 3;
   else imgftype = 2;

   attrib = MIN(atoi(argv[3]), NUMATTR-1);
   lambda = atof(argv[4]);
   if (argc>=6)  decision =  MIN(atoi(argv[5]), NUMDECISIONS-1);
   if (argc>=7)  outfname = argv[6];
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
   printf("Filtering image '%s' using attribute '%s' with lambda=%f\n", imgfname, Attribs[attrib].Name, lambda);
   printf("Decision rule: %s \n", Decisions[decision].Name);
   printf("Image: Width=%ld Height=%ld\n", width, height);
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
   
  NewAuxData = Attribs[attrib].NewData;
  DeleteAuxData = Attribs[attrib].DeleteData;
  AddToAuxData = Attribs[attrib].AddToData;  
  MergeAuxData = Attribs[attrib].MergeData;
  MergeToAuxData = Attribs[attrib].MergeToData;
  CloneAuxData = Attribs[attrib].CloneData;
   
  start = times(&tstruct);

  BuildTreeAndFilter();

  musec = (float)(times(&tstruct) - start)/((float)tickspersec);

  printf("wall-clock time: %f s\n",musec);
   
  if (imgftype==2) {
     r = ImagePGMBinWrite(outfname);
     free(out);  
  } else {
    if (VolStruct.DataType==1){
      for (i=0;i<size;i++)
        VolStruct.data.uchar_data[i]=out[i];
    } else {
      for (i=0;i<size;i++)
        VolStruct.data.ushort_data[i]=out[i];
    }
    free(out);  
    r = writeAVS(outfname, &avs_head, &VolStruct);
    DestroyVolume(&VolStruct);
  }
  if (r)  printf("Filtered image written to '%s'\n", outfname);

  for (i=0; i<size; i++) {
      if (node[i].Attribute) DeleteAuxData(node[i].Attribute);
  }
  free(remaining);
  free(reached);
  free(node);
  free(gval);
  return(0);
} /* main */
