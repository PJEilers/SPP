#include "TreeDiplomatic.h"

void RefTreeAreaFilterBerger(double lambda, greyval_t *out, greyval_t *gval, long size, Node *node_ref)
{
    
    pixel_t v, u, w, parent, lwb, upb;
    greyval_t val;
    
    lwb = 0;
    upb = size;
    for (v=lwb; v<upb; v++)
    {
        if (node_ref[v].isFiltered == false)
        {
	  w = v;
	  parent = node_ref[w].parent;
		
	  while ((parent != w) && (node_ref[w].isFiltered == false) &&
		 //((gval[w] == gval[parent]) || (node_ref[w].Area <= lambda)))
		 ((gval[w] == gval[parent]) || (node_ref[w].Area <= lambda) ))
            {
	      w = parent;
	      parent = node_ref[w].parent;
            }
	  if (node_ref[w].isFiltered == true)
            {
	      val = node_ref[w].filter;
	    }
	  else if (node_ref[w].Area > lambda)
            //else if (node_ref[w].attributes->area > lambda)
            {			
	      val = gval[w];
	    }
	  else
            {				
	      val = 0; // criterion cannot be satisfied
	    }
			
	  u = v;
	  while (u!=w)
            {
	      if ((lwb<=u) && (u<upb))
                {
		  node_ref[u].filter = val;
		  node_ref[u].isFiltered = true;
		}
	      u = node_ref[u].parent;
            }
	  if ((lwb<=w) && (w<upb))
            {
	      node_ref[w].filter = val;
	      node_ref[w].isFiltered = true;
	    }
        }
        
        out[v] = node_ref[v].filter;       
    }
}


void *runfilt(void *arg)
{	
	GeneralThreadData *thfiltdata = (GeneralThreadData *) arg;
        double lambda = thfiltdata->lambda;
    int self = thfiltdata->self;
    int nthreads = thfiltdata->nthreads;
    greyval_t *gval = thfiltdata->gval;
    greyval_t *out = thfiltdata->outRef;
    Node *node_ref= thfiltdata->node_ref;
    
    pixel_t v, u, w, parent, lwb, upb;
    greyval_t val;
    //lwb = 0;
    //upb = size;
    lwb = LWB(self, nthreads, thfiltdata->img.size2D, thfiltdata->img.depth);
    upb = UPB(self, nthreads, thfiltdata->img.size2D, thfiltdata->img.depth);
    
    printf("%d) lwb=%ld; upb=%ld.\n", self, lwb, upb);
    
    for (v=lwb; v<upb; v++)
    {
        if (node_ref[v].isFiltered == false)
        {
            w = v;
            parent = node_ref[w].parent;
			
			while ((parent != w) && (node_ref[w].isFiltered == false) &&
                    //((gval[w] == gval[parent]) || (node_ref[w].Area <= lambda)))
                    ((gval[w] == gval[parent]) || (node_ref[w].Area <= lambda) ))
            {
                w = parent;
                parent = node_ref[w].parent;
            }
            if (node_ref[w].isFiltered == true)
            {
				val = node_ref[w].filter;
			}
            //else if (node_ref[w].Area > lambda)
            else if (node_ref[w].Area > lambda)
            {			
				val = gval[w];
			}
            else
            {				
				val = 0; // criterion cannot be satisfied
			}
			
            u = v;
            while (u!=w)
            {
                if ((lwb<=u) && (u<upb))
                {
					node_ref[u].filter = val;
					node_ref[u].isFiltered = true;
				}
                u = node_ref[u].parent;
            }
            if ((lwb<=w) && (w<upb))
            {
				node_ref[w].filter = val;
				node_ref[w].isFiltered = true;
			}
        }
        
        out[v] = node_ref[v].filter;       
    }
    return NULL;
}

/*void FreeThreadFiltData(ThreadFiltData *data, int numthreads)
{
    free(data);
}*/

GeneralThreadData *MakeThreadFiltData(int numthreads, double lambda, GeneralThreadData *data)
{
    int i;
    for (i=0; i<numthreads; i++)
    {
        data[i].lambda = lambda;
    }
    return(data);
}

void RunFilter(GeneralThreadData *thdata, int nthreads)
{
    pthread_t threadID[MAXTHREADS];
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, runfilt, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

void ParallelFilter(int nthreads, GeneralThreadData *threadData)
{  
    RunFilter(threadData, threadData->numQTZLEVELS);
   // FreeThreadFiltData(thfiltdata, nthreadsRef);
}

void Filter(GeneralThreadData *threadData, double lambda)
{	
	RefTreeAreaFilterBerger(lambda, threadData->outRef, threadData->gval, threadData->img.size, threadData->node_ref);
	//ParallelFilter(nthreadsRef);
}
