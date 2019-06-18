#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <pthread.h>

#define false 0
#define true  1
#define MAXTHREADS 32
#define MAXLEVELS  256
#define CONNECTIVITY  4

typedef short bool;
typedef unsigned char ubyte;
/* typedef unsigned int uint; 
   typedef unsigned long ulong; */

#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))
#define LWB(self) ((self*size)/nthreads)
#define UPB(self) (((self+1)*size)/nthreads)
#define bottom (-1)

int nthreads;
pthread_t threadID[MAXTHREADS];

ubyte *gval=NULL, *out;
ulong width, height, size;

typedef struct MaxNode 
{ 
   ulong parent;
   ulong Area;
   void *Attribute;
   int filter; /* gray level after filtering */
} MaxNode;

MaxNode *node; 
bool *reached;

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

#include "decision.c"

/************************** queue ************************************/

typedef struct Queue
{
   ulong *Pixels;
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
   hq->Head = hq->Tail = 0;
   for (i=1; i<MAXLEVELS; i++)
   {
      hq[i].Pixels = hq[i-1].Pixels + numpixelsperlevel[i-1];
      hq[i].Head = hq[i].Tail = 0;
   }
   return(hq);
} /* QueueCreate */

#define QueueFirst(hq, h) (hq[h].Pixels[hq[h].Head++])
#define QueueAdd(hq, h, p) hq[h].Pixels[hq[h].Tail++] = p
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
      return(-2);
   }
   fscanf(infile, "P2\n");
   while ((c=fgetc(infile)) == '#')
      while ((c=fgetc(infile)) != '\n');
   ungetc(c, infile);
   fscanf(infile, "%lu %lu\n255\n", &width, &height);
   size = width*height;
   gval = malloc(size);
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
   int c;

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
   gval = malloc(size);
   if (gval==NULL) {
     fprintf (stderr, "Out of memory!");
     fclose(infile);
     return(-1);
   } 
   fread(gval, 1, size, infile);
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

   outfile = fopen(fname, "wb");
   if (outfile==NULL) {
      fprintf (stderr, "Error: Can't write the image: %s !", fname);
      return(-1);
   }
   fprintf(outfile, "P5\n%ld %ld\n255\n", width, height);
   fwrite(out, 1, (size_t)(size), outfile); 
 
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
   if ((p>=width) && (p-width>=lwb))    neighbors[n++] = p-width;
   if ((x>0) && (p-1>=lwb))             neighbors[n++] = p-1;
   if (p+width<upb)                     neighbors[n++] = p+width;
   return(n);
} /* GetNeighbors */

int BuildLocalTree(int self)
{
  ulong neighbors[CONNECTIVITY], numpixelsperlevel[MAXLEVELS], lero[MAXLEVELS];
  ulong p, q, x, y, lwb, upb, xm;
  int lev, fq;
  int numneighbors, m;
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
  QueueAdd(wait, lev, xm);
  reached[xm] = true;
  lero[lev] = xm;             

  while (lev >=0) {
    if (QueueNotEmpty(wait, lev)) { /* Encounter */
        node[lero[lev]].Area++;
        p = QueueFirst(wait, lev);
        if (gval[p]!=lev) printf("wrong!....");
        x = p % width; y = p / width;
        if (node[lero[lev]].Attribute) AddToAuxData(node[lero[lev]].Attribute, x, y); 
        else {
             node[lero[lev]].Attribute = NewAuxData(x, y);
             if (node[lero[lev]].Attribute==NULL)  return(-1);
	} 
        numneighbors = GetNeighbors(p, neighbors, lwb, upb);
        for (m=0; m<numneighbors; m++) {
            q = neighbors[m];
            if (!reached[q]) {
               fq = gval[q];            
               reached[q] = true;
               if (lero[fq]==bottom) lero[fq] = q;
               else node[q].parent = lero[fq];
               QueueAdd(wait, fq, q);
               if (fq > lev) lev = fq; 
            }
        }
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

ulong levroot(ulong x) {
  if (x==bottom) return bottom;
  while ((node[x].parent!=bottom) && (gval[x]==gval[node[x].parent]))
     x = node[x].parent;
  return x;
}

ulong Par(ulong x) {
  return levroot(node[x].parent);
}

void Connect (ulong x, ulong y) {

  ulong area = 0, area1 = 0;
  void *cor = NULL, *copa = NULL;
  ulong h, z;

  x = levroot(x);
  y = levroot(y);
  if (gval[y] > gval[x]) {
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

void levrootfix(ulong x, ulong lwb, ulong upb) {
  ulong h, u=x;
  x = levroot(x);
  while (u!=x && u>= lwb && u<upb) {
    h = node[u].parent;
    node[u].parent = x;
    u = h;
  }
}

void Fuse(int self, int i) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
  ulong lwb, upb;
  ulong p, x, y;
  int j;
   
  lwb = LWB(self);
  j = self+2*i-1; 
  upb = UPB(j); 
  if (upb>size) upb=size;

  /* get the horizontal boundary */
  j = self+i;
  x = LWB(j); 
  p = x % width;
  if ((p>0) && (x-1>=lwb)) {
     y = x-1;
     Connect(x, y);
  }  
  for (p=0; p<width && x<upb; p++){
      if ((x>=width) && (x-width>=lwb)) {
         y= x-width;
         Connect(x, y);
      }
      x++;
  }
  for (x=lwb; x<upb; x++)  levrootfix(x, lwb, upb);
}

/****************** Concurrent construction and filter of Maxtree  ****************************/

void *ccaf(void *arg) {
  int self= (int) arg, q, i;

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

/* checking the region being fused*/

  else {
    ulong neighbors[CONNECTIVITY];
    int numneighbors,j;
    ulong x,y,count=0;
    for (x=0; x<size; x++) {
        count++;
        if ((node[x].parent!=bottom) && (gval[x]<gval[node[x].parent])) {
           printf("something wrong in making parent!\n");
           break;
        }
        numneighbors = GetNeighbors(x, neighbors, 0, size);
        for (j=0; j<numneighbors; j++) {
            y = neighbors[j];
	    if ((gval[x]==gval[y]) && (node[x].parent!=node[y].parent) && !((node[x].parent==y) || (node[y].parent==x))) {
                printf("something wrong in fusing!\n");
                x = size;
                break;}
        }
    }
    printf("%lu pixels over %lu have been checked!\n",count,size);
  }


  Barrier(self);
  Decisions[decision].Filter(self, Attribs[attrib].Attribute, lambda);
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
   int r;
   ulong i;

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
   attrib = atoi(argv[3]);
   lambda = atof(argv[4]);
   if (argc>=6)  decision = atoi(argv[5]);
   if (argc>=7)  outfname = argv[6];
   
   ImagePGMRead(imgfname);
   if (gval==NULL) return(-1);
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
 
   out =  malloc(size);
   if (out==NULL) {
      fprintf(stderr, "Can't create output image! \n");
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
   
  BuildTreeAndFilter();
   
  r = ImagePGMBinWrite(outfname);
  if (!r)  printf("Filtered image written to '%s'\n", outfname);

  for (i=0; i<size; i++) {
      if (node[i].Attribute) DeleteAuxData(node[i].Attribute);
  }
  
  free(out);  
  free(reached);
  free(node);
  free(gval);
  return(0);
} /* main */
