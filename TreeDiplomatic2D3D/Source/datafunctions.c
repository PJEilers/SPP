#include "TreeDiplomatic.h"
void InitAuxDataStore( AuxDataStore *store, 
		       int itemsize, 
		       int maxchunks, 
		       int chunksize       )
{ 
  store->itemsize = itemsize;
  store->maxchunks = maxchunks;
  store->chunksize = chunksize;
  store->curchunk = 0;
  store->curitem = 0;
  store->data = calloc(maxchunks,sizeof(char*));
  store->data[0]= calloc(chunksize,itemsize);
}

void ClearAuxDataStore(AuxDataStore *store)
{ 
  int i;
  for ( i = 0; i<=store->curchunk; i++)
    free(store->data[i]);
  free(store->data);
}

void *GetNewAuxData(AuxDataStore *store)
{
  char *p;
  if (store->curitem <store->chunksize){
    p = store->data[store->curchunk] + store->itemsize * store->curitem;
  } else {
    store->curchunk ++;
    store->curitem = 0;
    store->data[store->curchunk] = 
      (char *)SafeCalloc(store->chunksize,store->itemsize);
    p = store->data[store->curchunk];
    if (p==NULL)
      printf("whoops\n");
  }
  store->curitem ++;

  return( (void *) p );    
}

/****** Typedefs and functions for area attributes ******************************/


void *NewAreaData( AuxDataStore *store, ulong x, ulong y, ulong z)
{
   AreaData *areadata;

   areadata = GetNewAuxData(store);
   areadata->Area = 1;
   return(areadata);
} /* NewAreaData */

void DeleteAreaData(void *areaattr)
{

} /* DeleteAreaData */

void AddToAreaData(void *areaattr, ulong x, ulong y, ulong z)
{
   AreaData *areadata = areaattr;

   areadata->Area ++;
} /* AddToAreaData */

void MergeAreaData(void *areaattr, void *childattr)
{
   AreaData *areadata = areaattr;
   AreaData *childdata = childattr;

   areadata->Area += childdata->Area;
} /* MergeAreaData */

void MergeToAreaData( AuxDataStore *store, void **thisattr, void *areaattr, void *childattr)
{
   AreaData *thisdata = *thisattr;
   AreaData *areadata = areaattr;
   AreaData *childdata = childattr;

   if (thisdata==NULL) {
      thisdata = GetNewAuxData(store);
      *thisattr = thisdata;
   }
   thisdata->Area = areadata->Area + childdata->Area;
} /* MergeToAreaData */

void CloneAreaData( AuxDataStore *store, void **thisattr, void *areaattr)
{
   AreaData *thisdata = *thisattr;
   AreaData *areadata = areaattr;

   if (thisdata==NULL) {
      thisdata = GetNewAuxData(store);
      *thisattr = thisdata;
   }
   thisdata->Area = areadata->Area;
} /* CloneAreaData */

double AreaAttribute(void *areaattr)
{
   AreaData *areadata = areaattr;
   double area;

   area = areadata->Area;
   return(area);
} /* AreaAttribute */

/****** Typedefs and functions for minimum enclosing rectangle attributes *******/

void *NewEnclRectData( AuxDataStore *store, ulong x, ulong y, ulong z)
{
   EnclRectData *rectdata;

   rectdata = GetNewAuxData(store);
   rectdata->MinX = rectdata->MaxX = x;
   rectdata->MinY = rectdata->MaxY = y;
   rectdata->MinZ = rectdata->MaxZ = z;
   return(rectdata);
} /* NewEnclRectData */

void DeleteEnclRectData(void *rectattr)
{
} /* DeleteEnclRectData */

void AddToEnclRectData(void *rectattr, ulong x, ulong y, ulong z)
{
   EnclRectData *rectdata = rectattr;

   rectdata->MinX = MIN(rectdata->MinX, x);
   rectdata->MinY = MIN(rectdata->MinY, y);
   rectdata->MinZ = MIN(rectdata->MinZ, z);
   rectdata->MaxX = MAX(rectdata->MaxX, x);
   rectdata->MaxY = MAX(rectdata->MaxY, y);
   rectdata->MaxZ = MAX(rectdata->MaxZ, z);
} /* AddToEnclRectData */

void MergeEnclRectData(void *rectattr, void *childattr)
{
   EnclRectData *rectdata = rectattr;
   EnclRectData *childdata = childattr;

   rectdata->MinX = MIN(rectdata->MinX, childdata->MinX);
   rectdata->MinY = MIN(rectdata->MinY, childdata->MinY);
   rectdata->MinZ = MIN(rectdata->MinZ, childdata->MinZ);
   rectdata->MaxX = MAX(rectdata->MaxX, childdata->MaxX);
   rectdata->MaxY = MAX(rectdata->MaxY, childdata->MaxY);
   rectdata->MaxZ = MAX(rectdata->MaxZ, childdata->MaxZ);
} /* MergeEnclRectData */

void MergeToEnclRectData( AuxDataStore *store, void **thisattr, void *rectattr, void *childattr)
{
   EnclRectData *thisdata = *thisattr;
   EnclRectData *rectdata = rectattr;
   EnclRectData *childdata = childattr;

   if (thisdata==NULL) {
      thisdata = GetNewAuxData(store);
      *thisattr = thisdata;
   }
   thisdata->MinX = MIN(rectdata->MinX, childdata->MinX);
   thisdata->MinY = MIN(rectdata->MinY, childdata->MinY);
   thisdata->MinZ = MIN(rectdata->MinZ, childdata->MinZ);
   thisdata->MaxX = MAX(rectdata->MaxX, childdata->MaxX);
   thisdata->MaxY = MAX(rectdata->MaxY, childdata->MaxY);
   thisdata->MaxZ = MAX(rectdata->MaxZ, childdata->MaxZ);
} /* MergeToEnclRectData */

void CloneEnclRectData( AuxDataStore *store, void **thisattr, void *rectattr)
{
   EnclRectData *thisdata = *thisattr;
   EnclRectData *rectdata = rectattr;

   if (thisdata==NULL) {
      thisdata = GetNewAuxData(store);
      *thisattr = thisdata;
   }
   thisdata->MinX = rectdata->MinX;
   thisdata->MinY = rectdata->MinY;
   thisdata->MinZ = rectdata->MinZ;
   thisdata->MaxX = rectdata->MaxX;
   thisdata->MaxY = rectdata->MaxY;
   thisdata->MaxZ = rectdata->MaxZ;
} /* CloneEnclRectData */

double EnclRectAreaAttribute(void *rectattr)
{
   EnclRectData *rectdata = rectattr;
   double volume;

   volume = (rectdata->MaxX - rectdata->MinX + 1)
          * (rectdata->MaxY - rectdata->MinY + 1)
          * (rectdata->MaxZ - rectdata->MinZ + 1);
   return(volume);
} /* EnclRectAreaAttribute */

double EnclRectDiagAttribute(void *rectattr)
/* Computes the square of the length of the diagonal */
{
   EnclRectData *rectdata = rectattr;
   double minx, miny, minz, maxx, maxy, maxz, l;

   minx = rectdata->MinX;
   miny = rectdata->MinY;
   minz = rectdata->MinZ;
   maxx = rectdata->MaxX;
   maxy = rectdata->MaxY;
   maxz = rectdata->MaxZ;
   l = (maxx-minx+1) * (maxx-minx+1)
     + (maxy-miny+1) * (maxy-miny+1)
     + (maxz-minz+1) * (maxz-minz+1);
   return(l);
} /* EnclRectDiagAttribute */

/****** Typedefs and functions for moment of inertia attributes **************************/


void *NewInertiaData( AuxDataStore *store, ulong x, ulong y, ulong z)  
{
   InertiaData *inertiadata;

   inertiadata = GetNewAuxData(store);
   inertiadata->Area = 1;
   inertiadata->SumX = x;
   inertiadata->SumY = y;
   inertiadata->SumZ = z;
   inertiadata->SumR2 = x*x + y*y + z*z;
   return(inertiadata);
} /* NewInertiaData */

void DeleteInertiaData(void *inertiaattr)
{
} /* DeleteInertiaData */

void AddToInertiaData(void *inertiaattr, ulong x, ulong y, ulong z)
{
   InertiaData *inertiadata = inertiaattr;

   inertiadata->Area ++;
   inertiadata->SumX += x;
   inertiadata->SumY += y;
   inertiadata->SumZ += z;
   inertiadata->SumR2 += x*x + y*y + z*z;
} /* AddToInertiaData */

void MergeInertiaData(void *inertiaattr, void *childattr)
{
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;

   inertiadata->Area += childdata->Area;
   inertiadata->SumX += childdata->SumX;
   inertiadata->SumY += childdata->SumY;
   inertiadata->SumZ += childdata->SumZ;
   inertiadata->SumR2 += childdata->SumR2;
} /* MergeInertiaData */

void MergeToInertiaData( AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr)
{
   InertiaData *thisdata = *thisattr;
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;

   if (thisdata==NULL) {
      thisdata = GetNewAuxData(store);
      *thisattr= thisdata;
   }
   thisdata->Area = inertiadata->Area + childdata->Area;
   thisdata->SumX = inertiadata->SumX + childdata->SumX;
   thisdata->SumY = inertiadata->SumY + childdata->SumY;
   thisdata->SumZ = inertiadata->SumZ + childdata->SumZ;
   thisdata->SumR2 = inertiadata->SumR2 + childdata->SumR2;
} /* MergeToInertiaData */

void CloneInertiaData( AuxDataStore *store, void **thisattr, void *inertiaattr)
{
   InertiaData *thisdata = *thisattr;
   InertiaData *inertiadata = inertiaattr;

   if (thisdata==NULL) {
      thisdata = GetNewAuxData(store);
      *thisattr = thisdata;
   }
   thisdata->Area = inertiadata->Area;
   thisdata->SumX = inertiadata->SumX;
   thisdata->SumY = inertiadata->SumY;
   thisdata->SumZ = inertiadata->SumZ;
   thisdata->SumR2 = inertiadata->SumR2;
} /* CloneInertiaData */

double InertiaAttribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double area, inertia;

   area = inertiadata->Area;
   inertia = inertiadata->SumR2 -
             (inertiadata->SumX * inertiadata->SumX +
              inertiadata->SumY * inertiadata->SumY +
              inertiadata->SumZ * inertiadata->SumZ) / area
             + area / 6.0;  /* ??? */
   return(inertia);
} /* InertiaAttribute */

double InertiaDivA2Attribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double inertia, area;

   area = (double)(inertiadata->Area);
   inertia = inertiadata->SumR2 -
             (inertiadata->SumX * inertiadata->SumX +
              inertiadata->SumY * inertiadata->SumY +
              inertiadata->SumZ * inertiadata->SumZ) / area
             + area / 6.0;  /* ??? */
   return(inertia/pow(area,5.0/3.0));
} /* InertiaDivA2Attribute */

double MeanXAttribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double area, sumx;

   area = inertiadata->Area;
   sumx = inertiadata->SumX;
   return(sumx/area);
} /* MeanXAttribute */

double MeanYAttribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double area, sumy;

   area = inertiadata->Area;
   sumy = inertiadata->SumY;
   return(sumy/area);
} /* MeanYAttribute */

double MeanZAttribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double area, sumz;

   area = inertiadata->Area;
   sumz = inertiadata->SumZ;
   return(sumz/area);
} /* MeanZAttribute */

static AttribStruct Attribs[NUMATTR] =
{
  {"Area", sizeof(AreaData), NewAreaData, DeleteAreaData, AddToAreaData, MergeAreaData, MergeToAreaData, CloneAreaData, AreaAttribute},
  {"Area of min. enclosing rectangle", sizeof(EnclRectData),  NewEnclRectData, DeleteEnclRectData, AddToEnclRectData, MergeEnclRectData, 
      MergeToEnclRectData, CloneEnclRectData, EnclRectAreaAttribute},
  {"Square of diagonal of min. enclosing rectangle", sizeof(EnclRectData), NewEnclRectData, DeleteEnclRectData, AddToEnclRectData, 
      MergeEnclRectData, MergeToEnclRectData, CloneEnclRectData, EnclRectDiagAttribute},
  {"Moment of Inertia", sizeof(InertiaData), NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MergeToInertiaData, 
      CloneInertiaData, InertiaAttribute},
  {"(Moment of Inertia) / (area)^2", sizeof(InertiaData), NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, 
      MergeToInertiaData, CloneInertiaData, InertiaDivA2Attribute},
  {"Mean X position", sizeof(InertiaData), NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MergeToInertiaData, 
      CloneInertiaData, MeanXAttribute},
  {"Mean Y position", sizeof(InertiaData), NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MergeToInertiaData, 
      CloneInertiaData, MeanYAttribute},
  {"Mean Z position", sizeof(InertiaData), NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MergeToInertiaData, 
      CloneInertiaData, MeanZAttribute}
};

AttribStruct* getAttribs() {
    return Attribs;
}

/* Free all thread data*/

void FreeThreadData(GeneralThreadData *threadData) {
    int i;
    
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
    for ( i = 0; i < threadData->nthreads; i++) {
        ClearAuxDataStore(&threadData[i].store);
        ClearAuxDataStore(&threadData[i].refstore);
    }

    
    free(threadData);
}


/* Initialize all thread data */

void InitializeThreadData(greyval_t *gval, pixel_t *SORTED, int nthreads, ImageProperties img, GeneralThreadData *threadData, double lambda, int attrib) {
    greyval_t *outRef =  malloc(img.size*sizeof(greyval_t));
    int i;
    AttribStruct *Attribs = getAttribs();
    for(i = 0; i < nthreads; i++) {
        threadData[i].SORTED = SORTED;
        threadData[i].gval = gval; 
        threadData[i].nthreads = nthreads;
        threadData[i].img = img;
        threadData[i].outRef = outRef;
        threadData[i].lambda = lambda;
        threadData[i].attrib = attrib;
        InitAuxDataStore(&(threadData[i].store), 
                         Attribs[attrib].size, 
                         img.size/img.size2D, img.size2D);
        InitAuxDataStore(&(threadData[i].refstore), 
                        Attribs[attrib].size, 
                        img.size/img.size2D, img.size2D);
        
    }
}

/*************************************************************************************************************************************/