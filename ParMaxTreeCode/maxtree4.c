/* maxtree4.c */
/* September 21, 2004  Erik R. Urbach */
/* Max-tree with a single scalar or vector attribute parameter, */
/* one or more reference images, and an optional template. */
/* Attribute: I/(A^2) (default) and */
/*            geometric/central/norm./invariant moments */
/* Decision: Subtractive (default) and Direct */
/* Similarity: Euclidean (default) and Difference */
/* Criterion: (d(tau(C),ref) >= epsilon) */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>



#define NUMLEVELS     256
#define CONNECTIVITY  4
#define PI 3.14159265358979323846



typedef short bool;
#define false 0
#define true  1

typedef unsigned char ubyte;
typedef unsigned int uint;
typedef unsigned long ulong;



typedef struct ImageGray ImageGray;

struct ImageGray
{
   ulong Width;
   ulong Height;
   ubyte *Pixmap;
};



typedef struct HQueue
{
   ulong *Pixels;
   ulong Head;
   ulong Tail; /* First free place in queue, or -1 when the queue is full */
} HQueue;



typedef struct MaxNode MaxNode;

struct MaxNode
{
   MaxNode *Parent;
   ulong *Children;
  int noChildren;
   ulong Area;
   void *Attribute;
   ubyte Level;
   ubyte NewLevel;  /* gray level after filtering */
};



/* Status stores the information of the pixel status: the pixel can be
 * NotAnalyzed, InTheQueue or assigned to node k at level h. In this
 * last case Status(p)=k. */
#define ST_NotAnalyzed  -1
#define ST_InTheQueue   -2

typedef struct MaxTree MaxTree;

struct MaxTree
{
   long *Status;
   ulong *NumPixelsBelowLevel;
   ulong *NumNodesAtLevel; /* Number of nodes C^k_h at level h */
   MaxNode **Nodes;
   void *(*NewAuxData)(ulong, ulong);
   void (*AddToAuxData)(void *, ulong, ulong);
   void (*MergeAuxData)(void *, void *);
   void (*DeleteAuxData)(void *);
};

int h_min;

void MaxTreeDelete(MaxTree *mt);



typedef struct InertiaData
{
   ulong Area;
   double SumX, SumY, SumX2, SumY2;
} InertiaData;

void *NewInertiaData(ulong x, ulong y)
{
   InertiaData *inertiadata;

   inertiadata = malloc(sizeof(InertiaData));
   if (inertiadata)
   {
      inertiadata->Area = 1;
      inertiadata->SumX = x;
      inertiadata->SumY = y;
      inertiadata->SumX2 = x*x;
      inertiadata->SumY2 = y*y;
   }
   return(inertiadata);
} /* NewInertiaData */

void DeleteInertiaData(void *inertiaattr)
{
   free(inertiaattr);
} /* DeleteInertiaData */

void AddToInertiaData(void *inertiaattr, ulong x, ulong y)
{
   InertiaData *inertiadata = inertiaattr;

   inertiadata->Area ++;
   inertiadata->SumX += x;
   inertiadata->SumY += y;
   inertiadata->SumX2 += x*x;
   inertiadata->SumY2 += y*y;
} /* AddToInertiaData */

void MergeInertiaData(void *inertiaattr, void *childattr)
{
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;

   inertiadata->Area += childdata->Area;
   inertiadata->SumX += childdata->SumX;
   inertiadata->SumY += childdata->SumY;
   inertiadata->SumX2 += childdata->SumX2;
   inertiadata->SumY2 += childdata->SumY2;
} /* MergeInertiaData */

void InertiaDivA2Attribute(void *inertiaattr, double *vec)
{
   InertiaData *inertiadata = inertiaattr;
   double inertia, area;

   area = (double)(inertiadata->Area);
   inertia = inertiadata->SumX2 + inertiadata->SumY2 -
             (inertiadata->SumX * inertiadata->SumX +
              inertiadata->SumY * inertiadata->SumY) / area
             + area / 6.0;
   *vec = inertia*2.0*PI/(area*area);
} /* InertiaDivA2Attribute */



void *NewMomentData(ulong x, ulong y)
{
   double *attr, *m;
   double h, v=1.0;
   int p, q;

   m = attr = calloc(10, sizeof(double));
   if (attr==NULL)  return(NULL);
   for (q=0; q<4; q++)
   {
      h = v;
      for (p=0; p<4-q; p++)
      {
         *(m++) = h;
         h *= x;
      }
      v *= y;
   }
   return(attr);
} /* NewMomentData */

void DeleteMomentData(void *m)
{
   free(m);
} /* DeleteMomentData */

void AddToMomentData(void *momentattr, ulong x, ulong y)
{
   double *m = momentattr;
   double h, v=1.0;
   int p, q;

   for (q=0; q<4; q++)
   {
      h = v;
      for (p=0; p<4-q; p++)
      {
         *(m++) += h;
         h *= x;
      }
      v *= y;
   }
} /* AddToMomentData */

void MergeMomentData(void *momentattr, void *childattr)
{
   double *momentdata = momentattr;
   double *childdata = childattr;
   int i;

   for (i=0; i<10; i++)  momentdata[i] += childdata[i];
} /* MergeMomentData */

void MomentGeometricAttribute(void *momentattr, double *vec)
{
   memcpy(vec, momentattr, 10*sizeof(double));
} /* MomentGeometricAttribute */

void MomentCentralAttribute(void *momentattr, double *vec)
{
   double *m = momentattr;
   double xbar, ybar;

   xbar = m[1]/m[0];
   ybar = m[4]/m[0];
   vec[0] = m[0];
   vec[1] = m[2] - xbar*m[1];
   vec[2] = m[7] - ybar*m[4];
   vec[3] = m[5] - ybar*m[1];
   vec[4] = m[3] - 3.0*xbar*m[2] + 2.0*m[1]*xbar*xbar;
   vec[5] = m[8] - 2.0*ybar*m[5] - xbar*m[7] + 2.0*ybar*ybar*m[1];
   vec[6] = m[6] - 2.0*xbar*m[5] - ybar*m[2] + 2.0*xbar*xbar*m[4];
   vec[7] = m[9] - 3.0*ybar*m[7] + 2.0*ybar*ybar*m[4];
} /* MomentCentralAttribute */

void MomentNormalizedAttribute(void *momentattr, double *vec)
{
   double *m = momentattr;
   double xbar, ybar, mu00o2, mu00o3, twoxbar2, twoybar2;

   xbar = m[1]/m[0];
   ybar = m[4]/m[0];
   mu00o2 = m[0]*m[0];
   mu00o3 = pow(m[0], 2.5);
   twoxbar2 = 2.0*xbar*xbar;
   twoybar2 = 2.0*ybar*ybar;
   vec[0] = (m[2] - xbar*m[1]) / mu00o2;
   vec[1] = (m[7] - ybar*m[4]) / mu00o2;
   vec[2] = (m[5] - ybar*m[1]) / mu00o2;
   vec[3] = (m[3] - 3.0*xbar*m[2] + twoxbar2*m[1]) / mu00o3;
   vec[4] = (m[8] - 2.0*ybar*m[5] - xbar*m[7] + twoybar2*m[1]) / mu00o3;
   vec[5] = (m[6] - 2.0*xbar*m[5] - ybar*m[2] + twoxbar2*m[4]) / mu00o3;
   vec[6] = (m[9] - 3.0*ybar*m[7] + twoybar2*m[4]) / mu00o3;
} /* MomentNormalizedAttribute */

void MomentInvariantsAttribute(void *momentattr, double *vec)
{
   double *m = momentattr;
   double xbar, ybar, mu00o2, mu00o3, twoxbar2, twoybar2;
   double nu02, nu11, nu20, nu03, nu12, nu21, nu30;

   xbar = m[1]/m[0];
   ybar = m[4]/m[0];
   mu00o2 = m[0]*m[0];
   mu00o3 = pow(m[0], 2.5);
   twoxbar2 = 2.0*xbar*xbar;
   twoybar2 = 2.0*ybar*ybar;

   nu20 = (m[2] - xbar*m[1]) / mu00o2;
   nu02 = (m[7] - ybar*m[4]) / mu00o2;
   nu11 = (m[5] - ybar*m[1]) / mu00o2;
   nu30 = (m[3] - 3.0*xbar*m[2] + twoxbar2*m[1]) / mu00o3;
   nu12 = (m[8] - 2.0*ybar*m[5] - xbar*m[7] + twoybar2*m[1]) / mu00o3;
   nu21 = (m[6] - 2.0*xbar*m[5] - ybar*m[2] + twoxbar2*m[4]) / mu00o3;
   nu03 = (m[9] - 3.0*ybar*m[7] + twoybar2*m[4]) / mu00o3;
   vec[0] = nu20 + nu02;
   vec[1] = pow((nu20 - nu02), 2.0) + 4.0*pow(nu11,2.0);
   vec[2] = pow((nu30 - 3.0*nu12), 2.0) + pow((3.0*nu21 - nu03), 2.0);
   vec[3] = pow((nu30 + nu12), 2.0) + pow((nu21 + nu03), 2.0);
   vec[4] = (nu30-3.0*nu12)*(nu30+nu12)*(pow((nu30+nu12),2.0) - 3.0*pow((nu21+nu03),2.0)) +
            (3.0*nu21-nu03)*(nu21+nu03)*(3.0*pow((nu30+nu12),2.0)-pow((nu21+nu03),2.0));
   vec[5] = (nu20-nu02) * (pow((nu30+nu12),2.0)-pow((nu21+nu03),2.0)) +
            4.0*nu11*(nu30+nu12)*(nu21+nu03);
   vec[6] = (3.0*nu21-nu03)*(nu30+nu12)*(pow((nu30+nu12),2.0) -
            3.0*pow((nu21+nu03),2.0)) +
            (3.0*nu12-nu30)*(nu21+nu03)*(3.0*pow((nu30+nu12),2.0)-pow((nu21+nu03),2.0));
} /* MomentInvariantsAttribute */

void MomentYapInvariantAttribute(void *momentattr, double *vec)
{
   double *m = momentattr;
   double xbar, ybar, mu00, mu20, mu02, mu11, mu30, mu12, mu21, mu03;
   double mu20min02, theta, a2, a3, b, b2, b3, c, c2, c3;

   xbar = m[1]/m[0];
   ybar = m[4]/m[0];
   mu00 = m[0];
   mu20 = m[2] - xbar*m[1];
   mu02 = m[7] - ybar*m[4];
   mu11 = m[5] - ybar*m[1];
   mu30 = m[3] - 3.0*xbar*m[2] + 2.0*m[1]*xbar*xbar;
   mu12 = m[8] - 2.0*ybar*m[5] - xbar*m[7] + 2.0*ybar*ybar*m[1];
   mu21 = m[6] - 2.0*xbar*m[5] - ybar*m[2] + 2.0*xbar*xbar*m[4];
   mu03 = m[9] - 3.0*ybar*m[7] + 2.0*ybar*ybar*m[4];
   mu20min02 = mu20 - mu02;
   if (mu11==0.0)  theta = (mu20min02 >= 0.0) ? 0.0 : (-0.5*PI);
   else if (mu20min02==0.0)  theta = (mu11<0.0) ? (-0.25*PI) : (0.25*PI);
   else
   {
      theta = 0.5*atan(2.0*mu11/mu20min02);
      if (mu20min02<0.0)  theta += (mu11>0.0) ? (0.5*PI) : (-0.5*PI);
   }
   /*printf("theta = %f  mu30=%f mu03=%f\n", theta*180.0/PI, mu30, mu03);*/
   a2 = 1.0/(m[0]*m[0]);
   a3 = pow(m[0],-2.5);
   b = cos(theta);
   b2 = b*b;
   b3 = b2*b;
   c = sin(theta);
   c2 = c*c;
   c3 = c2*c;
   vec[0] = a2 * (b2*mu20 + 2.0*b*c*mu11 + c2*mu02);
   vec[1] = a2 * (-b*c*mu20 + (b2-c2)*mu11 + b*c*mu02);
   vec[2] = a2 * (c2*mu20 - 2.0*b*c*mu11 + b2*mu02);
   vec[3] = a3 * (b3*mu30 + 3.0*b2*c*mu21 + 3.0*b*c2*mu12 + c3*mu03);
   vec[4] = a3 * (-b2*c*mu30 + (b3-2.0*b*c2)*mu21 + (2.0*b2*c-c3)*mu12 + b*c2*mu03);
   vec[5] = a3 * (b*c2*mu30 + (c3-2.0*b2*c)*mu21 + (b3-2.0*b*c2)*mu12 + b2*c*mu03);
   vec[6] = a3 * (-c3*mu30 + 3.0*b*c2*mu21 - 3.0*b2*c*mu12 + b3*mu03);
} /* MomentYapInvariantAttribute */



int BinCoef[16] = {1, 1, 1, 1, 1, 2, 3, 1, 3, 1};

double Bin(int n, int k)
{
   return(BinCoef[k*4+n]);
} /* Bin */

void MomentYapTildeInvariantAttribute(void *momentattr, double *vec)
{
   double nu[16];
   double *m = momentattr;
   double xbar, ybar, mu00, mu20, mu02, mu11, mu30, mu12, mu21, mu03;
   double mu20min02, theta, a2, a3, b, b2, b3, c, c2, c3;
   double N=1000.0;
   int mm,n,p,q;

   xbar = m[1]/m[0];
   ybar = m[4]/m[0];
   mu00 = m[0];
   mu20 = m[2] - xbar*m[1];
   mu02 = m[7] - ybar*m[4];
   mu11 = m[5] - ybar*m[1];
   mu30 = m[3] - 3.0*xbar*m[2] + 2.0*m[1]*xbar*xbar;
   mu12 = m[8] - 2.0*ybar*m[5] - xbar*m[7] + 2.0*ybar*ybar*m[1];
   mu21 = m[6] - 2.0*xbar*m[5] - ybar*m[2] + 2.0*xbar*xbar*m[4];
   mu03 = m[9] - 3.0*ybar*m[7] + 2.0*ybar*ybar*m[4];
   mu20min02 = mu20 - mu02;
   if (mu11==0.0)  theta = (mu20min02 >= 0.0) ? 0.0 : (-0.5*PI);
   else if (mu20min02==0.0)  theta = (mu11<0.0) ? (-0.25*PI) : (0.25*PI);
   else
   {
      theta = 0.5*atan(2.0*mu11/mu20min02);
      if (mu20min02<0.0)  theta += (mu11>0.0) ? (0.5*PI) : (-0.5*PI);
   }
   printf("theta = %f  mu30=%f mu03=%f\n", theta*180.0/PI, mu30, mu03);
   a2 = 1.0/(m[0]*m[0]);
   a3 = pow(m[0],-2.5);
   b = cos(theta);
   b2 = b*b;
   b3 = b2*b;
   c = sin(theta);
   c2 = c*c;
   c3 = c2*c;
   nu[0] = 1.0;
   nu[1] = nu[4] = 0.0;
   nu[2] = a2 * (b2*mu20 + 2.0*b*c*mu11 + c2*mu02);
   nu[3] = a2 * (-b*c*mu20 + (b2-c2)*mu11 + b*c*mu02);
   nu[5] = a2 * (c2*mu20 - 2.0*b*c*mu11 + b2*mu02);
   nu[6] = a3 * (b3*mu30 + 3.0*b2*c*mu21 + 3.0*b*c2*mu12 + c3*mu03);
   nu[8] = a3 * (-b2*c*mu30 + (b3-2.0*b*c2)*mu21 + (2.0*b2*c-c3)*mu12 + b*c2*mu03);
   nu[9] = a3 * (b*c2*mu30 + (c3-2.0*b2*c)*mu21 + (b3-2.0*b*c2)*mu12 + b2*c*mu03);
   nu[12] = a3 * (-c3*mu30 + 3.0*b*c2*mu21 - 3.0*b2*c*mu12 + b3*mu03);
   for (mm=0; mm<4; mm++)
   {
      for (n=0; n<(4-mm); n++)
      {
         *vec = 0.0;
         for (p=0; p<=n; p++)
         {
            for (q=0; q<=mm; q++)
            {
               *vec += Bin(n,p)*Bin(mm,q)*pow(0.5*N*N,((p+q)/2.0)+1.0)*pow(0.5*N,n+mm-p-q)*nu[q*4+p];
            }
         }
         vec++;
      }
   }
} /* MomentYapTildeInvariantAttribute */



int Fac[4] = {1, 1, 2, 6};

double Poch(int a, int k)
{
  /* indexing */
  int i;

  /* initialisation */
  double result = 1.0;

  if (k<0)
    {
      printf("Error in Pochhammer calculation: (%d)%d \n",a,k);
      exit(0);
    }
  
  for (i=0; i<k; i++)
    {
      result *= (a+i); 
    }

  return result;
}

double rho(int n, double p, long N)
{
   return(pow(-1.0,n) * pow((1.0-p)/p,n) * (Fac[n]/Poch(-N,n)));
}

double coefficient(int k, int n, double p, int N)
{
  double coeff = 0.0;

  if (k<=0) 
    {
      coeff = 1.0;
    }
  else if (k>n) 
    {
      coeff = 0.0;
    }
  else if (n==1 && k==1)
    {
      coeff = -1.0/(N*p);
    }
  else  
    {
      coeff = coefficient(k,n-1,p,N) + (((n-1.0)*(1.0-p))/(p*(N-(n-1.0)))*(coefficient(k,n-1,p,N) - coefficient(k,n-2,p,N)) - (coefficient(k-1,n-1,p,N)/(p*(N-(n-1.0)))));
    }
  return coeff;
}

void MomentKrawtchoukAttribute(void *momentattr, double *vec)
{
   double nu[16], tilde[16];
   double *m = momentattr;
   double xbar, ybar, mu00, mu20, mu02, mu11, mu30, mu12, mu21, mu03;
   double mu20min02, theta, a2, a3, b, b2, b3, c, c2, c3;
   double p1=0.5, p2=0.5;
   double N = 10000.0;
   int mm,n,p,q,i,j;

   xbar = m[1]/m[0];
   ybar = m[4]/m[0];
   mu00 = m[0];
   mu20 = m[2] - xbar*m[1];
   mu02 = m[7] - ybar*m[4];
   mu11 = m[5] - ybar*m[1];
   mu30 = m[3] - 3.0*xbar*m[2] + 2.0*m[1]*xbar*xbar;
   mu12 = m[8] - 2.0*ybar*m[5] - xbar*m[7] + 2.0*ybar*ybar*m[1];
   mu21 = m[6] - 2.0*xbar*m[5] - ybar*m[2] + 2.0*xbar*xbar*m[4];
   mu03 = m[9] - 3.0*ybar*m[7] + 2.0*ybar*ybar*m[4];
   mu20min02 = mu20 - mu02;
   if (mu11==0.0)  theta = (mu20min02 >= 0.0) ? 0.0 : (-0.5*PI);
   else if (mu20min02==0.0)  theta = (mu11<0.0) ? (-0.25*PI) : (0.25*PI);
   else
   {
      theta = 0.5*atan(2.0*mu11/mu20min02);
      if (mu20min02<0.0)  theta += (mu11>0.0) ? (0.5*PI) : (-0.5*PI);
   }
   a2 = 1.0/(m[0]*m[0]);
   a3 = pow(m[0],-2.5);
   b = cos(theta);
   b2 = b*b;
   b3 = b2*b;
   c = sin(theta);
   c2 = c*c;
   c3 = c2*c;
   nu[0] = 1.0;
   nu[1] = nu[4] = 0.0;
   nu[2] = a2 * (b2*mu20 + 2.0*b*c*mu11 + c2*mu02);
   nu[3] = a2 * (-b*c*mu20 + (b2-c2)*mu11 + b*c*mu02);
   nu[5] = a2 * (c2*mu20 - 2.0*b*c*mu11 + b2*mu02);
   nu[6] = a3 * (b3*mu30 + 3.0*b2*c*mu21 + 3.0*b*c2*mu12 + c3*mu03);
   nu[8] = a3 * (-b2*c*mu30 + (b3-2.0*b*c2)*mu21 + (2.0*b2*c-c3)*mu12 + b*c2*mu03);
   nu[9] = a3 * (b*c2*mu30 + (c3-2.0*b2*c)*mu21 + (b3-2.0*b*c2)*mu12 + b2*c*mu03);
   nu[12] = a3 * (-c3*mu30 + 3.0*b*c2*mu21 - 3.0*b2*c*mu12 + b3*mu03);
   for (mm=0; mm<4; mm++)
   {
      for (n=0; n<(4-mm); n++)
      {
         tilde[mm*4+n] = 0.0;
         for (p=0; p<=n; p++)
         {
            for (q=0; q<=mm; q++)
            {
               tilde[mm*4+n] += Bin(n,p)*Bin(mm,q)*pow(0.5*N*N,((p+q)/2.0)+1.0)*pow(0.5*N,n+mm-p-q)*nu[q*4+p];
            }
         }
      }
   }
   for (mm=0; mm<4; mm++)
   {
      for (n=0; n<(4-mm); n++)
      {
         *vec = 0.0;
         for (i=0; i<=n; i++)
         {
            for (j=0; j<=mm; j++)
            {
               *vec += coefficient(i,n,p1,N)*coefficient(j,mm,p2,N)*tilde[j*4+i];
            }
         }
         *vec /= sqrt(rho(n,p1,N)*rho(mm,p2,N));
         vec++;
      }
   }
} /* MomentKrawtchoukAttribute */



double SimilarityEuclidean(double *obj, double *ref, int veclen)
{
   double d, sum=0.0;
   int i;

   for (i=0; i<veclen; i++)
   {
      d = obj[i] - ref[i];
      sum += d*d;
   }
   printf("dist=%f\n", sqrt(sum));
   return(sqrt(sum));
} /* SimilarityEuclidean */



double SimilarityDifference(double *obj, double *ref, int veclen)
{
   double sum=0.0;
   int i;

   for (i=0; i<veclen; i++)  sum += obj[i] - ref[i];
   return(sum);
} /* SimilarityDifference */



void PrintVector(double *vec, int veclen)
{
   int i;

   for (i=0; i<veclen; i++)  printf("%f ", vec[i]);
   printf("\n");
} /* PrintVector */



bool CriterionGEq(double *obj, int numref, double **ref, int veclen,
                  int k, double epsilon,
                  double (*similarity)(double *, double *, int))
{
   int i, num=0;

   for (i=0; i<numref; i++)
   {
      if ((*similarity)(obj, ref[i], veclen) < epsilon)
      {
         num++;
         if (num>=k)  return(false);
      }
   }
   return(true);
} /* CriterionGEq */



ImageGray *ImageGrayCreate(ulong width, ulong height)
{
   ImageGray *img;

   img = malloc(sizeof(ImageGray));
   if (img==NULL)  return(NULL);
   img->Width = width;
   img->Height = height;
   img->Pixmap = malloc(width*height);
   if (img->Pixmap==NULL)
   {
      free(img);
      return(NULL);
   }
   return(img);
} /* ImageGrayCreate */



void ImageGrayDelete(ImageGray *img)
{
   free(img->Pixmap);
   free(img);
} /* ImageGrayDelete */



void ImageGrayInit(ImageGray *img, ubyte h)
{
   memset(img->Pixmap, h, (img->Width)*(img->Height));
} /* ImageGrayInit */



ImageGray *ImagePGMBinRead(char *fname)
{
   FILE *infile;
   ImageGray *img;
   ulong width, height;
   int c;

   infile = fopen(fname, "rb");
   if (infile==NULL)  return(NULL);
   fscanf(infile, "P5\n");
   while ((c=fgetc(infile)) == '#')
      while ((c=fgetc(infile)) != '\n');
   ungetc(c, infile);
   fscanf(infile, "%lu %lu\n255\n", &width, &height);
   img = ImageGrayCreate(width, height);
   if (img)  fread(img->Pixmap, 1, width*height, infile);
   fclose(infile);
   return(img);
} /* ImagePGMBinRead */



int ImagePGMBinWrite(ImageGray *img, char *fname)
{
   FILE *outfile;

   outfile = fopen(fname, "wb");
   if (outfile==NULL)  return(-1);
   fprintf(outfile, "P5\n%ld %ld\n255\n", img->Width, img->Height);
   fwrite(img->Pixmap, 1, (size_t)((img->Width)*(img->Height)), outfile);
   fclose(outfile);
   return(0);
} /* ImagePGMBinWrite */



HQueue *HQueueCreate(ulong imgsize, ulong *numpixelsperlevel)
{
   HQueue *hq;
   int i;

   hq = calloc(NUMLEVELS, sizeof(HQueue));
   if (hq==NULL)  return(NULL);
   hq->Pixels = calloc(imgsize, sizeof(ulong));
   if (hq->Pixels==NULL)
   {
      free(hq);
      return(NULL);
   }
   hq->Head = hq->Tail = 0;
   for (i=1; i<NUMLEVELS; i++)
   {
      hq[i].Pixels = hq[i-1].Pixels + numpixelsperlevel[i-1];
      hq[i].Head = hq[i].Tail = 0;
   }
   return(hq);
} /* HQueueCreate */



void HQueueDelete(HQueue *hq)
{
   free(hq->Pixels);
   free(hq);
} /* HQueueDelete */



#define HQueueFirst(hq,h)  (hq[h].Pixels[hq[h].Head++])
#define HQueueAdd(hq,h,p)  hq[h].Pixels[hq[h].Tail++] = p
#define HQueueNotEmpty(hq,h)  (hq[h].Head != hq[h].Tail)



int GetNeighbors(ubyte *shape, ulong imgwidth, ulong imgsize, ulong p,
                 ulong *neighbors)
{
   ulong x;
   int n=0;

   x = p % imgwidth;
   if ((x<(imgwidth-1)) && (shape[p+1]))      neighbors[n++] = p+1;
   if ((p>=imgwidth) && (shape[p-imgwidth]))  neighbors[n++] = p-imgwidth;
   if ((x>0) && (shape[p-1]))                 neighbors[n++] = p-1;
   p += imgwidth;
   if ((p<imgsize) && (shape[p]))             neighbors[n++] = p;
   return(n);
} /* GetNeighbors */

void ArrayCopy(MaxNode *input, MaxNode *output, int size)
{
  int i;
  for (i=0;i<size;i++)
    {
      output[i] = input[i];
    }
} /* ArrayCopy */

void CopyIndices(ulong *input, ulong *output, int size)
{
  int i;
  for (i=0; i<size; i++)
    {
      output[i] = input[i];
    }

} /* CopyIndices */

void AddChild(MaxNode *node, int idx)
{
  ulong *result;
  int i;
  if (node->noChildren==NULL)
    {
      //      fprintf(stderr,"NULL Children\n");
      result = calloc (1,sizeof(ulong));
      //fprintf(stderr,"Memory assigned\n");
      result[0] = idx;
      //fprintf(stderr,"result[0] set\n");
      node->noChildren = 1;
      //fprintf(stderr,"node->noChildren set\n");
      node->Children = result;
      //fprintf(stderr,"Result copied\n");
    } else {
      //fprintf(stderr,"Existing Children\n");
      result = calloc(node->noChildren +1,sizeof(ulong));
      //fprintf(stderr,"Additional memory assigned\n");
      CopyIndices(node->Children,result,node->noChildren);
      //fprintf(stderr,"Indices copied\n");
      node->noChildren = node->noChildren +1;
      //fprintf(stderr,"node->noChildren increased t: %i\n",node->noChildren+1);
      free(node->Children);
      //fprintf(stderr,"Olde children array freed\n");
      node->Children = calloc(node->noChildren +1,sizeof(ulong));
      CopyIndices(result,node->Children,node->noChildren);
      free(result);
      //fprintf(stderr,"New children array copied\n");
    }
} /* AddChild */


int MaxTreeFlood(MaxTree *mt, HQueue *hq, ulong *numpixelsperlevel,
                 bool *nodeatlevel, ImageGray *img, ubyte *shape, int h,
                 ulong *thisarea,
                 void **thisattr)
/* Returns value >=NUMLEVELS if error */
{
   ulong neighbors[CONNECTIVITY];
   ubyte *pixmap;
   void *attr = NULL, *childattr;
   ulong imgwidth, imgsize, p, q, idx, x, y;
   ulong area = *thisarea, childarea;
   MaxNode *node,*parNode;
   int numneighbors, i;
   int m;

   imgwidth = img->Width;
   imgsize = imgwidth * (img->Height);
   pixmap = img->Pixmap;
   while(HQueueNotEmpty(hq, h))
   {
      area++;
      p = HQueueFirst(hq, h);
      x = p % imgwidth;
      y = p / imgwidth;
      if (attr)  mt->AddToAuxData(attr, x, y);
      else
      {
         attr = mt->NewAuxData(x, y);
         if (attr==NULL)  return(NUMLEVELS);
         if (*thisattr)  mt->MergeAuxData(attr, *thisattr);
      }
      mt->Status[p] = mt->NumNodesAtLevel[h];
      numneighbors = GetNeighbors(shape, imgwidth, imgsize, p, neighbors);
      for (i=0; i<numneighbors; i++)
      {
         q = neighbors[i];
         if (mt->Status[q]==ST_NotAnalyzed)
         {
            HQueueAdd(hq, pixmap[q], q);
            mt->Status[q] = ST_InTheQueue;
            nodeatlevel[pixmap[q]] = true;
            if (pixmap[q] > pixmap[p])
            {
               m = pixmap[q];
               childarea = 0;
               childattr = NULL;
               do
               {  
                  fprintf(stderr,"Value m: %i, h: %i\n",m,h);
                  fprintf(stderr,"Neighbour q: %i\n",q);

                  m = MaxTreeFlood(mt,hq,numpixelsperlevel,nodeatlevel,img,shape,m, &childarea, &childattr);
                  fprintf(stderr,"Value m: %i\n",m);
                  if (m>=NUMLEVELS)
                  {
                     mt->DeleteAuxData(attr);
                     return(m);
                  }
               } while (m!=h);
               area += childarea;
               mt->MergeAuxData(attr, childattr);
            }
         }
      }
   }
   fprintf(stderr,"Adding Node at level %i\n",h);
   node = (mt->Nodes[h]) + mt->NumNodesAtLevel[h];
   mt->NumNodesAtLevel[h] = mt->NumNodesAtLevel[h]+1;
   m = h-1;
   while ((m>=0) && (nodeatlevel[m]==false))  m--;
   if (m>=0)
   {
     fprintf(stderr,"Parent clause\n");
     //realloc(mt->Nodes[h],mt->NumNodesAtLevel[h]*sizeof(MaxNode));
     //fprintf(stderr,"Reallocated Memory\n");
     fprintf(stderr,"Assigned *node\n");
     node->Parent = (mt->Nodes[m]) + mt->NumNodesAtLevel[m]-1;
     fprintf(stderr,"End parent clause\n");

   } else {
     fprintf(stderr,"No parent clause\n");
     //realloc(mt->Nodes[h],mt->NumNodesAtLevel[h]*sizeof(MaxNode));
     node->Parent = node;
     fprintf(stderr,"End no parent clause\n");
   }
   //fprintf(stderr,"Assigning attributes\n");
   node->Area = area;
   //fprintf(stderr,"Area\n");
   node->Attribute = attr;
   //fprintf(stderr,"Attributes\n");
   node->Level = h;
   //fprintf(stderr,"Level\n");
   nodeatlevel[h] = false;
   //fprintf(stderr,"nodeatlevel\n");
   *thisarea = area;
   //fprintf(stderr,"*thisarea\n");
   *thisattr = attr;
   //fprintf(stderr,"*thisattr\n");
   //fprintf(stderr,"Value m: %i\n",m);
   return(m);
} /* MaxTreeFlood */



MaxTree *MaxTreeCreate(ImageGray *img, ImageGray *template,
                       void *(*newauxdata)(ulong, ulong),
                       void (*addtoauxdata)(void *, ulong, ulong),
                       void (*mergeauxdata)(void *, void *),
                       void (*deleteauxdata)(void *))
{
   ulong numpixelsperlevel[NUMLEVELS];
   bool nodeatlevel[NUMLEVELS];
   HQueue *hq;
   MaxTree *mt;
   ubyte *pixmap = img->Pixmap;
   void *attr = NULL;
   ulong imgsize, p, m=0, area=0;
   int l,i;

   /* Allocate structures */
   mt = malloc(sizeof(MaxTree));
   if (mt==NULL)  return(NULL);
   imgsize = (img->Width)*(img->Height);
   mt->Status = calloc((size_t)imgsize, sizeof(long));
   if (mt->Status==NULL)
   {
      free(mt);
      return(NULL);
   }
   mt->NumPixelsBelowLevel = calloc(NUMLEVELS, sizeof(ulong));
   if (mt->NumPixelsBelowLevel==NULL)
   {
      free(mt->Status);
      free(mt);
      return(NULL);
   }
   mt->NumNodesAtLevel = calloc(NUMLEVELS, sizeof(ulong));
   if (mt->NumNodesAtLevel==NULL)
   {
      free(mt->NumPixelsBelowLevel);
      free(mt->Status);
      free(mt);
      return(NULL);
   }
   mt->Nodes = calloc((size_t)NUMLEVELS, sizeof(MaxNode*));
   for (i=0;i<NUMLEVELS;i++)
     {
       mt->Nodes[i] = calloc(imgsize/2 ,sizeof(MaxNode));
     }
   if (mt->Nodes==NULL)
   {
      free(mt->NumNodesAtLevel);
      free(mt->NumPixelsBelowLevel);
      free(mt->Status);
      free(mt);
      return(NULL);
   }

   /* Initialize structures */
   for (p=0; p<imgsize; p++)  mt->Status[p] = ST_NotAnalyzed;
   bzero(nodeatlevel, NUMLEVELS*sizeof(bool));
   bzero(numpixelsperlevel, NUMLEVELS*sizeof(ulong));
   /* Following bzero is redundant, array is initialized by calloc */
   /* bzero(mt->NumNodesAtLevel, NUMLEVELS*sizeof(ulong)); */
   for (p=0; p<imgsize; p++)  numpixelsperlevel[pixmap[p]]++;
   mt->NumPixelsBelowLevel[0] = 0;
   for (l=1; l<NUMLEVELS; l++)
   {
      mt->NumPixelsBelowLevel[l] = mt->NumPixelsBelowLevel[l-1] + numpixelsperlevel[l-1];
   }
   hq = HQueueCreate(imgsize, numpixelsperlevel);
   if (hq==NULL)
   {
      free(mt->Nodes);
      free(mt->NumNodesAtLevel);
      free(mt->NumPixelsBelowLevel);
      free(mt->Status);
      free(mt);
      return(NULL);
   }

   /* Find pixel m which has the lowest intensity l in the image */
   for (p=0; p<imgsize; p++)
   {
      if (pixmap[p]<pixmap[m])  m = p;
   }
   l = pixmap[m];
   h_min = l;

   /* Add pixel m to the queue */
   nodeatlevel[l] = true;
   HQueueAdd(hq, l, m);
   mt->Status[m] = ST_InTheQueue;

   /* Build the Max-tree using a flood-fill algorithm */
   mt->NewAuxData = newauxdata;
   mt->AddToAuxData = addtoauxdata;
   mt->MergeAuxData = mergeauxdata;
   mt->DeleteAuxData = deleteauxdata;
   l = MaxTreeFlood(mt, hq, numpixelsperlevel, nodeatlevel, img,
                    template->Pixmap, l, &area, &attr);
   fprintf(stderr,"MaxTreeFlooded\n");
   if (l>=NUMLEVELS)  MaxTreeDelete(mt);
   HQueueDelete(hq);
   return(mt);
} /* MaxTreeCreate */



void MaxTreeDelete(MaxTree *mt)
{
   void *attr;
   ulong i;
   int h;

   for (h=0; h<NUMLEVELS; h++)
   {
      for (i=0; i<mt->NumNodesAtLevel[h]; i++)
      {
         attr = mt->Nodes[h][i].Attribute;
         if (attr)  mt->DeleteAuxData(attr);
      }
   }
   free(mt->Nodes);
   free(mt->NumNodesAtLevel);
   free(mt->NumPixelsBelowLevel);
   free(mt->Status);
   free(mt);
} /* MaxTreeDelete */

/* int MaxTreeFilterDirect(MaxTree *mt, ImageGray *img, ImageGray *template,
                        ImageGray *out, void (*attribute)(void *, double *),
                        double (*similarity)(double *, double *, int),
                        bool criterion(double *, int, double **, int, int,
                                double,
                                double (*similarity)(double *, double *, int)),
                        int veclen, int numref, double **ref, int k,
                        double epsilon)
{
   MaxNode *node;
   double *vec;
   ubyte *shape = template->Pixmap;
   ulong i, idx, parent;
   int l;

   vec = calloc(veclen, sizeof(double));
   if (vec==NULL)  return(-1);
   for (l=0; l<NUMLEVELS; l++)
   {
      for (i=0; i<mt->NumNodesAtLevel[l]; i++)
      {
         idx = mt->NumPixelsBelowLevel[l] + i;
         node = &(mt->Nodes[idx]);
         parent = node->Parent;
         (*attribute)(node->Attribute, vec);
         if ((idx==parent) || ((*criterion)(vec,numref,ref,veclen,k,epsilon,similarity)))
         {
            node->NewLevel = node->Level;
         } else  node->NewLevel = mt->Nodes[parent].NewLevel;
      }
   }
   for (i=0; i<(img->Width)*(img->Height); i++)
   {
      if (shape[i])
      {
         idx = mt->NumPixelsBelowLevel[img->Pixmap[i]] + mt->Status[i];
         out->Pixmap[i] = mt->Nodes[idx].NewLevel;
      }
   }
   free(vec);
   return(0);
} 

*/
/*
int MaxTreeFilterSubtractive(MaxTree *mt, ImageGray *img, ImageGray *template,
                           ImageGray *out, void (*attribute)(void *, double *),
                           double (*similarity)(double *, double *, int),
                           bool criterion(double *, int, double **, int, int,
                                double,
                                double (*similarity)(double *, double *, int)),
                           int veclen, int numref, double **ref, int k,
                           double epsilon)
{
   MaxNode *node, *parnode;
   double *vec;
   ubyte *shape = template->Pixmap;
   ulong i, idx, parent;
   int l;

   vec = calloc(veclen, sizeof(double));
   if (vec==NULL)  return(-1);
   for (l=0; l<NUMLEVELS; l++)
   {
      for (i=0; i<mt->NumNodesAtLevel[l]; i++)
      {
         idx = mt->NumPixelsBelowLevel[l] + i;
         node = &(mt->Nodes[idx]);
         parent = node->Parent;
         parnode = &(mt->Nodes[parent]);
         (*attribute)(node->Attribute, vec);
         printf("l=%d i=%ld   tau(C) : ", l, i);
         PrintVector(vec, veclen);
         if ((idx==parent) || ((*criterion)(vec,numref,ref,veclen,k,epsilon,similarity)))
         {
            node->NewLevel = ((int)(node->Level)) + ((int)(parnode->NewLevel)) - ((int)(parnode->Level));
         } else  node->NewLevel = parnode->NewLevel;
      }
   }
   for (i=0; i<(img->Width)*(img->Height); i++)
   {
      if (shape[i])
      {
         idx = mt->NumPixelsBelowLevel[img->Pixmap[i]] + mt->Status[i];
         out->Pixmap[i] = mt->Nodes[idx].NewLevel;
      }
   }
   free(vec);
   return(0);
   } 
*/



ImageGray *GetTemplate(char *templatefname, ImageGray *img)
{
   ImageGray *template;

   if (templatefname)
   {
      template = ImagePGMBinRead(templatefname);
      if (template==NULL)  return(NULL);
      if ((img->Width != template->Width) || (img->Height != template->Height))
      {
	 ImageGrayDelete(template);
         return(NULL);
      }
   } else {
      template = ImageGrayCreate(img->Width, img->Height);
      if (template)  ImageGrayInit(template, NUMLEVELS-1);
   }
   return(template);
} /* GetTemplate */



int ComputeReferenceVector(char *reffname,
                           void (*attribute)(void *, double *),
                           double *ref)
{
   ImageGray *img;
   void *attr=NULL;
   ulong x, y;

   img = ImagePGMBinRead(reffname);
   if (img==NULL)  return(-1);
   for (y=0; y<(img->Height); y++)
   {
      for (x=0; x<(img->Width); x++)
      {
         if (img->Pixmap[y*(img->Width)+x])
         {
            if (attr)  AddToMomentData(attr, x, y);
            else
            {
               attr = NewMomentData(x, y);
               if (attr==NULL)
               {
                  ImageGrayDelete(img);
                  return(-1);
	       }
	    }
         }
      }
   }
   attribute(attr, ref);
   DeleteMomentData(attr);
   ImageGrayDelete(img);
   return(0);
} /* ComputeReferenceVector */



double **CreateReferenceVectors(int numref, char **reffname, int veclen,
                                void (*attribute)(void *, double *))
{
   double **ref;
   int i, r;

   ref = calloc(numref, sizeof(double *));
   if (ref==NULL)  return(NULL);
   ref[0] = calloc(veclen*numref, sizeof(double));
   if (ref[0]==NULL)
   {
      free(ref);
      return(NULL);
   }
   for (i=1; i<numref; i++)  ref[i] = ref[i-1]+veclen;
   for (i=0; i<numref; i++)
   {
      r = ComputeReferenceVector(reffname[i], attribute, ref[i]);
      if (r)
      {
         free(ref[0]);
         free(ref);
         return(NULL);
      }
   }
   return(ref);
} /* CreateReferenceVectors */



void DeleteReferenceVectors(double **ref)
{
   free(ref[0]);
   free(ref);
} /* DeleteReferenceVectors */



typedef struct AttrTp
{
   void (*Attribute)(void *, double *);
   char *Name;
} AttrTp;

int NumAttributes = 7;
AttrTp Attributes[7] =
{
   {MomentGeometricAttribute, "Geometric moments"},
   {MomentCentralAttribute, "Central moments"},
   {MomentNormalizedAttribute, "Normalized central moments"},
   {MomentInvariantsAttribute, "Hu's moment invariants"},
   {MomentYapInvariantAttribute, "Moment invariants (Yap)"},
   {MomentYapTildeInvariantAttribute, "Modified Moment Invariants (Yap)"},
   {MomentKrawtchoukAttribute, "Krawtchouk moment invariants (Yap)"}
};

int main(int argc, char *argv[])
{
   ImageGray *img, *template, *out;
   MaxTree *mt;
   char *imgfname, *templatefname = NULL, **reffnames, *outfname = "out.pgm";
   double **ref, epsilon=0.01;
   int i,r,veclen=10,numref,k=1, attrib=3;

   if (argc<3)
   {
      printf("Usage: %s <input image> <numref> {references} [k] [epsilon] [attrib] [output image] [template]\n", argv[0]);
      printf("Where attrib is:\n");
      for (i=0; i<NumAttributes; i++)  printf("\t%d : %s\n", i, Attributes[i].Name);
      exit(0);
   }
   imgfname = argv[1];
   numref = atoi(argv[2]);
   reffnames = argv+3;
   if (argc>=(4+numref))  k = atoi(argv[3+numref]);
   if (argc>=(5+numref))  epsilon = atof(argv[4+numref]);
   if (argc>=(6+numref))  attrib = atoi(argv[5+numref]);
   if (argc>=(7+numref))  outfname = argv[6+numref];
   if (argc>=(8+numref))  templatefname = argv[7+numref];
   printf("Using %s attribute\n", Attributes[attrib].Name);
   ref = CreateReferenceVectors(numref, reffnames, veclen, Attributes[attrib].Attribute);
   if (ref==NULL)
   {
      fprintf(stderr, "Can't read reference images\n");
      return(-1);
   }
   printf("Reference vectors:\n");
   for (i=0; i<numref; i++)
   {
      printf("%d : ", i);
      PrintVector(ref[i], veclen);
   }
   img = ImagePGMBinRead(imgfname);
   if (img==NULL)
   {
      fprintf(stderr, "Can't read image '%s'\n", imgfname);
      return(-1);
   }
   template = GetTemplate(templatefname, img);
   if (template==NULL)
   {
      fprintf(stderr, "Can't create template\n");
      ImageGrayDelete(img);
      return(-1);
   }
   out = ImageGrayCreate(img->Width, img->Height);
   if (out==NULL)
   {
      fprintf(stderr, "Can't create output image\n");
      ImageGrayDelete(template);
      ImageGrayDelete(img);
      return(-1);
   }
   mt = MaxTreeCreate(img, template, NewMomentData, AddToMomentData, MergeMomentData, DeleteMomentData);
   fprintf(stderr,"MaxTree created\n");
   if (mt==NULL)
   {
      fprintf(stderr, "Can't create Max-tree\n");
      ImageGrayDelete(out);
      ImageGrayDelete(template);
      ImageGrayDelete(img);
      return(-1);
   }
   //r = MaxTreeFilterDirect(mt, img, template, out, Attributes[attrib].Attribute, SimilarityEuclidean, CriterionGEq, veclen, numref, ref, k, epsilon);
   MaxTreeDelete(mt);
   r = ImagePGMBinWrite(out, outfname);
   if (r)  fprintf(stderr, "Error writing image '%s'\n", outfname);
   ImageGrayDelete(out);
   ImageGrayDelete(template);
   ImageGrayDelete(img);
   DeleteReferenceVectors(ref);
   return(0);
} /* main */
