#include<math.h>

typedef struct{
  int itemsize, 
    chunksize, 
    maxchunks,
    curchunk, 
    curitem;
  char **data;
} AuxDataStore;

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

typedef struct AreaData
{
   ulong Area;
} AreaData;

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

typedef struct EnclRectData 
{
   ulong MinX;
   ulong MinY;
   ulong MinZ;
   ulong MaxX;
   ulong MaxY;
   ulong MaxZ;
} EnclRectData;

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

typedef struct InertiaData 
{
   ulong Area;
   double SumX, SumY, SumZ, SumR2;
} InertiaData;

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

void MaxTreeFilterMin(int self, double (*attribute)(void *), double lambda, long *fdepth)
{
   ulong v, u, w, r, parent, lwb, upb;
   int val;
   
   lwb = LWB(self); upb = UPB(self);
   for (v=lwb; v<upb; v++) {
      if (node[v].filter == -1) {
         w = levroot(v,fdepth);
         r = w;
         parent = Par(w, fdepth);
         while ((parent != bottom) && (node[w].filter == -1)){
               if ((*attribute)(node[w].Attribute) < lambda) r=w;
               w = parent;
               parent = Par(w,fdepth);
         }
         if ((node[w].filter != -1) && (node[w].filter != gval[w])) { val = node[w].filter; r=w; }
         else if ((*attribute)(node[r].Attribute) >= lambda) val = gval[r]; 
         else if (((*attribute)(node[r].Attribute) < lambda) && (node[r].parent != bottom)) val = gval[node[r].parent]; 
         else val = 0; /* criterion cannot be satisfied */
         u = v;
         while (u!=r) {
           if ((lwb<=u) && (u<upb)) node[u].filter = val;
           u = node[u].parent;
           }
         if ((lwb<=r) && (r<upb)) node[r].filter = val;
       }
       out[v] = node[v].filter;
    }
} /* MaxTreeFilterMin */

void MaxTreeFilterDirect(int self, double (*attribute)(void *), double lambda, long *fdepth)
{
   ulong v, u, w, parent, lwb, upb;
   int val;
 
   lwb = LWB(self); upb = UPB(self);
   for (v=lwb; v<upb; v++) {
      if (node[v].filter == -1) {
         w = v;
         parent = node[w].parent;
         while ((parent != bottom) && (node[w].filter == -1) && 
            ((gval[w] == gval[parent]) || ((*attribute)(node[w].Attribute) < lambda))) {
                  w = parent;
                  parent = node[w].parent;
         }
         if (node[w].filter != -1) val = node[w].filter;
         else if ((*attribute)(node[w].Attribute) >= lambda) val = gval[w]; 
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
} /* MaxTreeFilterDirect */


void MaxTreeFilterMax(int self, double (*attribute)(void *), double lambda, long *fdepth) 
{
   ulong v, u, w, parent, lwb, upb;
   int val;
 
   lwb = LWB(self); upb = UPB(self);
   for (v=lwb; v<upb; v++) {
      if (node[v].filter == -1) {
         w = v;
         parent = node[w].parent;
         while ((parent != bottom) && (node[w].filter == -1) && 
            ((gval[w] == gval[parent]) || ((*attribute)(node[w].Attribute) < lambda))) {
                  w = parent;
                  parent = node[w].parent;
         }
         if (node[w].filter != -1) val = node[w].filter;
         else if ((*attribute)(node[w].Attribute) >= lambda) {
              val = gval[w]; 
              u = parent;
              while ((u != bottom) && (node[u].filter!=gval[u])) {
                 if ((lwb<=u) && (u<upb)) node[u].filter = gval[u];
                 else {remaining[u]=true; break;}      /* else let thread 0 do the remainning work */
                 u = node[u].parent;
	      }
         }
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
} /* MaxTreeFilterMax */

void MaxTreeFilterSubtractive(int self, double (*attribute)(void *), double lambda, long *fdepth) /* ????*/
{
   ulong v, u, w, lwb, upb;
   ulong levpath[MAXLEVELS];
   int i, j, filter[MAXLEVELS];

   lwb = LWB(self); upb = UPB(self);
   for (v=lwb; v<upb; v++) {
      if (node[v].filter == -1) {
         w=levroot(v,fdepth);       
         i=0;
         while ((w != bottom) && (node[w].filter == -1)) {
            levpath[i++] = w; 
            w = Par(w,fdepth);
         } 
         levpath[i] = w;
         if (w !=bottom) filter[i] = node[w].filter;
         for (j=i-1; j>=0; j--) {
             u = levpath[j];
             if (w!=bottom) {
                if (((*attribute)(node[u].Attribute)) < lambda) filter[j] = filter[j+1];
                else filter[j] = gval[u]+filter[j+1]-gval[w];
             } else if (((*attribute)(node[u].Attribute)) >= lambda)
                filter[j] = gval[u];
	     else filter[j]=0; 
	     w = u;
         }
         u = v;
         for (j=0; j<=i; j++) {
	   while (u!=levpath[j]) {
              if ((lwb<=u) && (u<upb)) node[u].filter = filter[j];
              u = node[u].parent;
           }
           if ((lwb<=u) && (u<upb)) node[u].filter = filter[j];
           if (u!=bottom) u = node[u].parent;
         }
      }
       out[v] = node[v].filter;
   }
} /* MaxTreeFilterSubtractive */

void *(*NewAuxData)( AuxDataStore *, ulong, ulong, ulong);
void (*DeleteAuxData)(void *);
void (*AddToAuxData)(void *, ulong, ulong, ulong);
void (*MergeAuxData)(void *, void *);
void (*MergeToAuxData)( AuxDataStore *, void **, void *, void *);
void (*CloneAuxData)( AuxDataStore *, void **, void *);

typedef struct AttribStruct
{
   char *Name;
   int size;
   void *(*NewData)( AuxDataStore *, ulong, ulong, ulong);
   void (*DeleteData)(void *);
   void (*AddToData)(void *, ulong, ulong, ulong);
   void (*MergeData)(void *, void *);
   void (*MergeToData)( AuxDataStore *, void **, void *, void *);
   void (*CloneData)( AuxDataStore *, void **, void *);
   double (*Attribute)(void *);
} AttribStruct;

#define NUMATTR 8

AttribStruct Attribs[NUMATTR] =
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

typedef struct DecisionStruct
{
   char *Name;
   void (*Filter)(int, double (*attribute)(void *), double, long *);

} DecisionStruct;

#define NUMDECISIONS 4
DecisionStruct Decisions[NUMDECISIONS] =
{
   {"Min", MaxTreeFilterMin},
   {"Direct", MaxTreeFilterDirect},
   {"Max", MaxTreeFilterMax},
   {"Subtractive", MaxTreeFilterSubtractive},
};

