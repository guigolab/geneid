/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: SortSR                                             *
*                                                              *
*   Sort by position start the SR set                          *
***************************************************************/

#include "blast2gff.h"

int Select2SR(sr_t* x, sr_t* y)
{
  int a;

  if (x == NULL)
    a = 1;
  else
    if (y == NULL)
      a = 0;
    else
      if (x->Start < y->Start)
	a = 0;
      else
	a = 1;
  
  return(a);
}

int SelectSR(sr_t* c[])
{
  int x,y,z,r;

  x = Select2SR(c[0],c[1]);
  
  y = Select2SR(c[2],c[3]);
  y = y + 2;

  z = Select2SR(c[4],c[5]);
  z = z + 4;

  /* All the possible options are... */
  if (c[x] == NULL && c[y] == NULL && c[z] == NULL)
    r = NOTFOUND;

  if (c[x] != NULL && c[y] == NULL && c[z] == NULL)
    r = x;

  if (c[x] == NULL && c[y] != NULL && c[z] == NULL)
    r = y;

  if (c[x] == NULL && c[y] == NULL && c[z] != NULL)
    r = z;

  if (c[x] != NULL && c[y] != NULL && c[z] == NULL)
    r = (c[x]->Start < c[y]->Start)? x : y;

  if (c[x] != NULL && c[y] == NULL && c[z] != NULL)
    r = (c[x]->Start < c[z]->Start)? x : z;

  if (c[x] == NULL && c[y] != NULL && c[z] != NULL)
    r = (c[y]->Start < c[z]->Start)? y : z;

  if (c[x] != NULL && c[y] != NULL && c[z] != NULL)
    {
      r = (c[x]->Start < c[y]->Start)? x : y;
      r = (c[z]->Start < c[r]->Start)? z : r;
    }

  return(r);
} 

void SortSR(packSR *allSR)
{
  int i;
  sr_t* c[STRANDS*FRAMES];
  long n;

  /* merge-sorting all the strands and frames */ 
  c[0] = allSR->sr[0];
  c[1] = allSR->sr[1];
  c[2] = allSR->sr[2];
  c[3] = allSR->sr[3];
  c[4] = allSR->sr[4];
  c[5] = allSR->sr[5];

  n = 0;
  /* i=0..5 and i=NOTFOUND means finish of the process */
  i = SelectSR(c);
  
  while (i != NOTFOUND)
    {
      /* Sharing pointers... */
      allSR->sortSR[n] = c[i]; 
      n++;

      /* Next SR of this list */
      c[i] = c[i]->next;
      i = SelectSR(c);
    }

  allSR->nSR = n;
}

