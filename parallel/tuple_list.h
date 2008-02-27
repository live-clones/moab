/*------------------------------------------------------------------------------
  
  Tuple list definition and utilities
  
  Conceptually, a tuple list is a list of n records or tuples,
  each with mi integers, ml longs, and mr reals
  (these types are defined in "types.h" as sint, slong, real;
   it may be that sint==slong)
  
  There are three arrays, one for each type (vi,vl,vr),
  with records layed out contiguously within each array

  ----------------------------------------------------------------------------*/

#ifndef TUPLE_LIST_H
#define TUPLE_LIST_H

/* requires "errmem.h" and "types.h" */
#if !defined(ERRMEM_H) || !defined(TYPES_H) || !defined(MINMAX_H) || !defined(SORT_H)
#warning "tuple_list.h" requires "errmem.h" and "types.h" and  "minmax.h" and "sort.h"
#endif

typedef struct {
  unsigned mi,ml,mul,mr;
  uint n, max;
  sint *vi; slong *vl; ulong *vul; real *vr;
} tuple_list;

/* storage layed out as: vi[max][mi], vl[max][ml], vr[max][mr]
   where "tuple" i is given by (vi[i][0:mi-1],vl[i][0:ml-1],vr[i][0:mr-1]).
   only the first n tuples are in use */

static void tuple_list_init_max(tuple_list *tl,
                                unsigned mi, unsigned ml, unsigned mul, 
                                unsigned mr, uint max)
{
  tl->n=0; tl->max=max;
  tl->mi=mi,tl->ml=ml,tl->mul=mul,tl->mr=mr;
  tl->vi=tmalloc(sint, max*mi);
  tl->vl=tmalloc(slong,max*ml);
  tl->vul=tmalloc(ulong,max*mul);
  tl->vr=tmalloc(real, max*mr);
}

static void tuple_list_free(tuple_list *tl) {
  free(tl->vi), free(tl->vl), free(tl->vul), free(tl->vr);
}

static void tuple_list_resize(tuple_list *tl, uint max)
{
  tl->max = max;
  tl->vi=trealloc(sint, tl->vi,tl->max*tl->mi);
  tl->vl=trealloc(slong,tl->vl,tl->max*tl->ml);
  tl->vul=trealloc(ulong,tl->vul,tl->max*tl->mul);
  tl->vr=trealloc(real, tl->vr,tl->max*tl->mr);
}

static void tuple_list_grow(tuple_list *tl)
{
  tuple_list_resize(tl,tl->max+tl->max/2+1);
}

void tuple_list_permute(tuple_list *tl, uint *perm, void *work);
void tuple_list_sort(tuple_list *tl, unsigned key, buffer *buf);

#endif

