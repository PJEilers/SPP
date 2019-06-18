#define PI 3.14159265358979323846
int attrib, decision=1; /*  the original default value of decision is 3; */
double lambda;

/****** Typedefs and functions for area attributes ******************************/

typedef struct AreaData
{
   ulong Area;
} AreaData;

void *NewAreaData(ulong x, ulong y)
{
   AreaData *areadata;

   areadata = SafeMalloc(sizeof(AreaData));
   areadata->Area = 1;
   return(areadata);
} /* NewAreaData */

void DeleteAreaData(void *areaattr)
{
   free(areaattr);
} /* DeleteAreaData */

void AddToAreaData(void *areaattr, ulong x, ulong y)
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

void MergeToAreaData(void **thisattr, void *areaattr, void *childattr)
{
   AreaData *thisdata = *thisattr;
   AreaData *areadata = areaattr;
   AreaData *childdata = childattr;

   if (thisdata==NULL) {
      thisdata = SafeMalloc(sizeof(AreaData));
      *thisattr = thisdata;
   }
   thisdata->Area = areadata->Area + childdata->Area;
} /* MergeToAreaData */

void CloneAreaData(void **thisattr, void *areaattr)
{
   AreaData *thisdata = *thisattr;
   AreaData *areadata = areaattr;

   if (thisdata==NULL) {
      thisdata = SafeMalloc(sizeof(AreaData));
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
   ulong MaxX;
   ulong MaxY;
} EnclRectData;

void *NewEnclRectData(ulong x, ulong y)
{
   EnclRectData *rectdata;

   rectdata = SafeMalloc(sizeof(EnclRectData));
   rectdata->MinX = rectdata->MaxX = x;
   rectdata->MinY = rectdata->MaxY = y;
   return(rectdata);
} /* NewEnclRectData */

void DeleteEnclRectData(void *rectattr)
{
   free(rectattr);
} /* DeleteEnclRectData */

void AddToEnclRectData(void *rectattr, ulong x, ulong y)
{
   EnclRectData *rectdata = rectattr;

   rectdata->MinX = MIN(rectdata->MinX, x);
   rectdata->MinY = MIN(rectdata->MinY, y);
   rectdata->MaxX = MAX(rectdata->MaxX, x);
   rectdata->MaxY = MAX(rectdata->MaxY, y);
} /* AddToEnclRectData */

void MergeEnclRectData(void *rectattr, void *childattr)
{
   EnclRectData *rectdata = rectattr;
   EnclRectData *childdata = childattr;

   rectdata->MinX = MIN(rectdata->MinX, childdata->MinX);
   rectdata->MinY = MIN(rectdata->MinY, childdata->MinY);
   rectdata->MaxX = MAX(rectdata->MaxX, childdata->MaxX);
   rectdata->MaxY = MAX(rectdata->MaxY, childdata->MaxY);
} /* MergeEnclRectData */

void MergeToEnclRectData(void **thisattr, void *rectattr, void *childattr)
{
   EnclRectData *thisdata = *thisattr;
   EnclRectData *rectdata = rectattr;
   EnclRectData *childdata = childattr;

   if (thisdata==NULL) {
      thisdata = SafeMalloc(sizeof(EnclRectData));
      *thisattr = thisdata;
   }
   thisdata->MinX = MIN(rectdata->MinX, childdata->MinX);
   thisdata->MinY = MIN(rectdata->MinY, childdata->MinY);
   thisdata->MaxX = MAX(rectdata->MaxX, childdata->MaxX);
   thisdata->MaxY = MAX(rectdata->MaxY, childdata->MaxY);
} /* MergeToEnclRectData */

void CloneEnclRectData(void **thisattr, void *rectattr)
{
   EnclRectData *thisdata = *thisattr;
   EnclRectData *rectdata = rectattr;

   if (thisdata==NULL) {
      thisdata = SafeMalloc(sizeof(EnclRectData));
      *thisattr = thisdata;
   }
   thisdata->MinX = rectdata->MinX;
   thisdata->MinY = rectdata->MinY;
   thisdata->MaxX = rectdata->MaxX;
   thisdata->MaxY = rectdata->MaxY;
} /* CloneEnclRectData */

double EnclRectAreaAttribute(void *rectattr)
{
   EnclRectData *rectdata = rectattr;
   double area;

   area = (rectdata->MaxX - rectdata->MinX + 1)
        * (rectdata->MaxY - rectdata->MinY + 1);
   return(area);
} /* EnclRectAreaAttribute */

double EnclRectDiagAttribute(void *rectattr)
/* Computes the square of the length of the diagonal */
{
   EnclRectData *rectdata = rectattr;
   double minx, miny, maxx, maxy, l;

   minx = rectdata->MinX;
   miny = rectdata->MinY;
   maxx = rectdata->MaxX;
   maxy = rectdata->MaxY;
   l = (maxx-minx+1) * (maxx-minx+1)
     + (maxy-miny+1) * (maxy-miny+1);
   return(l);
} /* EnclRectDiagAttribute */



/****** Typedefs and functions for moment of inertia attributes **************************/

typedef struct InertiaData 
{
   ulong Area;
   double SumX, SumY, SumX2, SumY2;
} InertiaData;

void *NewInertiaData(ulong x, ulong y)
{
   InertiaData *inertiadata;

   inertiadata = SafeMalloc(sizeof(InertiaData));
   inertiadata->Area = 1;
   inertiadata->SumX = x;
   inertiadata->SumY = y;
   inertiadata->SumX2 = x*x;
   inertiadata->SumY2 = y*y;
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

void MergeToInertiaData(void **thisattr, void *inertiaattr, void *childattr)
{
   InertiaData *thisdata = *thisattr;
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;

   if (thisdata==NULL) {
      thisdata = SafeMalloc(sizeof(InertiaData));
      *thisattr= thisdata;
   }
   thisdata->Area = inertiadata->Area + childdata->Area;
   thisdata->SumX = inertiadata->SumX + childdata->SumX;
   thisdata->SumY = inertiadata->SumY + childdata->SumY;
   thisdata->SumX2 = inertiadata->SumX2 + childdata->SumX2;
   thisdata->SumY2 = inertiadata->SumY2 + childdata->SumY2;
} /* MergeToInertiaData */

void CloneInertiaData(void **thisattr, void *inertiaattr)
{
   InertiaData *thisdata = *thisattr;
   InertiaData *inertiadata = inertiaattr;

   if (thisdata==NULL) {
      thisdata = SafeMalloc(sizeof(InertiaData));
      *thisattr = thisdata;
   }
   thisdata->Area = inertiadata->Area;
   thisdata->SumX = inertiadata->SumX;
   thisdata->SumY = inertiadata->SumY;
   thisdata->SumX2 = inertiadata->SumX2;
   thisdata->SumY2 = inertiadata->SumY2;
} /* CloneInertiaData */

double InertiaAttribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double area, inertia;

   area = inertiadata->Area;
   inertia = inertiadata->SumX2 + inertiadata->SumY2 -
             (inertiadata->SumX * inertiadata->SumX +
              inertiadata->SumY * inertiadata->SumY) / area
             + area / 6.0;
   return(inertia);
} /* InertiaAttribute */

double InertiaDivA2Attribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double inertia, area;

   area = (double)(inertiadata->Area);
   inertia = inertiadata->SumX2 + inertiadata->SumY2 -
             (inertiadata->SumX * inertiadata->SumX +
              inertiadata->SumY * inertiadata->SumY) / area
             + area / 6.0;
   return(inertia*2.0*PI/(area*area));
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

void MaxTreeFilterMin(int self, double (*attribute)(void *), double lambda)
{
   ulong v, u, w, r, parent;
   int val;

   for (v=LWB(self); v<UPB(self); v++) {
      if (node[v].filter == -1) {
         w = v;
         parent = node[w].parent;
         if ((parent != bottom) && (gval[w] == gval[parent])) {
            w = parent;
            parent = node[w].parent;
         }
         r = w;
         while ((parent != bottom) && (node[w].filter == -1)){
               if ((*attribute)(node[w].Attribute) < lambda) r=w;
               w = parent;
               parent = node[w].parent;
         }
         if ((node[w].filter != -1) && (node[w].filter != gval[w])) { val = node[w].filter; r=w; }
         else if ((*attribute)(node[r].Attribute) >= lambda) val = gval[r]; 
         else if (((*attribute)(node[r].Attribute) < lambda) && (node[r].parent != bottom)) val = gval[node[r].parent]; 
         else val = 0; /* criterion cannot be satisfied */
         u = v;
         while (u!=r) {
           if ((LWB(self)<=u) && (u<UPB(self))) node[u].filter = val;
           u = node[u].parent;
           }
         if ((LWB(self)<=r) && (r<UPB(self))) node[r].filter = val;
       }
       out[v] = node[v].filter;
    }
} /* MaxTreeFilterMin */

void MaxTreeFilterDirect(int self, double (*attribute)(void *), double lambda)
{
   ulong v, u, w, parent;
   int val;

   for (v=LWB(self); v<UPB(self); v++) {
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
           if ((LWB(self)<=u) && (u<UPB(self))) node[u].filter = val;
           u = node[u].parent;
           }
         if ((LWB(self)<=w) && (w<UPB(self))) node[w].filter = val;
       }
       out[v] = node[v].filter;
   }
} /* MaxTreeFilterDirect */


void MaxTreeFilterMax(int self, double (*attribute)(void *), double lambda) /* ????*/
{
   ulong v, u, w, parent;
   int val;

   for (v=LWB(self); v<UPB(self); v++) {
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
           if ((LWB(self)<=u) && (u<UPB(self))) node[u].filter = val;
           u = node[u].parent;
           }
         if ((LWB(self)<=w) && (w<UPB(self))) node[w].filter = val;
       }
       out[v] = node[v].filter;
   }
} /* MaxTreeFilterMax */

void MaxTreeFilterSubtractive(int self, double (*attribute)(void *), double lambda) /* ????*/
{
   ulong v, u, w, r, parent;
   int val;

   for (v=LWB(self); v<UPB(self); v++) {
      if (node[v].filter == -1) {
         w = v;
         parent = node[w].parent;
         if ((parent != bottom) && (gval[w] == gval[parent])) {
            w = parent;
            parent = node[w].parent;
         }
         r = w;
         while ((parent != bottom) && (node[w].filter == -1)){
               if ((*attribute)(node[w].Attribute) < lambda) r=w;
               w = parent;
               parent = node[w].parent;
         }
         if ((node[w].filter != -1) && (node[w].filter != gval[w])) { val = node[w].filter; r=w; }
         else if ((*attribute)(node[r].Attribute) >= lambda) val = gval[r]; 
         else if (((*attribute)(node[r].Attribute) < lambda) && (node[r].parent != bottom)) val = gval[node[r].parent]; 
         else val = 0; /* criterion cannot be satisfied */
         u = v;
         while (u!=r) {
           if ((LWB(self)<=u) && (u<UPB(self))) node[u].filter = val;
           u = node[u].parent;
           }
         if ((LWB(self)<=r) && (r<UPB(self))) node[r].filter = val;
       }
       out[v] = node[v].filter;
    }
} /* MaxTreeFilterSubtractive */

void *(*NewAuxData)(ulong, ulong);
void (*DeleteAuxData)(void *);
void (*AddToAuxData)(void *, ulong, ulong);
void (*MergeAuxData)(void *, void *);
void (*MergeToAuxData)(void **, void *, void *);
void (*CloneAuxData)(void **, void *);

typedef struct AttribStruct
{
   char *Name;
   void *(*NewData)(ulong, ulong);
   void (*DeleteData)(void *);
   void (*AddToData)(void *, ulong, ulong);
   void (*MergeData)(void *, void *);
   void (*MergeToData)(void **, void *, void *);
   void (*CloneData)(void **, void *);
   double (*Attribute)(void *);
} AttribStruct;

#define NUMATTR 7

AttribStruct Attribs[NUMATTR] =
{
  {"Area", NewAreaData, DeleteAreaData, AddToAreaData, MergeAreaData, MergeToAreaData, CloneAreaData, AreaAttribute},
  {"Area of min. enclosing rectangle", NewEnclRectData, DeleteEnclRectData, AddToEnclRectData, MergeEnclRectData, 
      MergeToEnclRectData, CloneEnclRectData, EnclRectAreaAttribute},
  {"Square of diagonal of min. enclosing rectangle", NewEnclRectData, DeleteEnclRectData, AddToEnclRectData, 
      MergeEnclRectData, MergeToEnclRectData, CloneEnclRectData, EnclRectDiagAttribute},
  {"Moment of Inertia", NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MergeToInertiaData, 
      CloneInertiaData, InertiaAttribute},
  {"(Moment of Inertia) / (area)^2", NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, 
      MergeToInertiaData, CloneInertiaData, InertiaDivA2Attribute},
  {"Mean X position", NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MergeToInertiaData, 
      CloneInertiaData, MeanXAttribute},
  {"Mean Y position", NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MergeToInertiaData, 
      CloneInertiaData, MeanYAttribute}
};

typedef struct DecisionStruct
{
   char *Name;
   void (*Filter)(int, double (*attribute)(void *), double);

} DecisionStruct;

#define NUMDECISIONS 4
DecisionStruct Decisions[NUMDECISIONS] =
{
   {"Min", MaxTreeFilterMin},
   {"Direct", MaxTreeFilterDirect},
   {"Max", MaxTreeFilterMax},
   {"Subtractive", MaxTreeFilterSubtractive},
};

