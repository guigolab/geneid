/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: SortHSP                                            *
*                                                              *
*   Sort by position startQ the HSP set                        *
***************************************************************/

#include "blast2gff.h"

int Split(hsp** hsps,int i, int j)
{
  long k;
  int x,c1,c2,m;
  int pivot;
  int pivotStartQ;
  hsp* tmp;
  
  /* Choosing the random pivot */
  k = hsps[i]->StartQ;

  /* How many elements have less StartQ value than pivot */
  for(m = 0, x=i+1; x <=j; x++)
    if (hsps[x]->StartQ <= k)
	m++;
  
  /* This is the right place for the pivot */
  pivot = m+i;

  tmp = hsps[i];
  hsps[i] = hsps[pivot];
  hsps[pivot] = tmp;

  pivotStartQ = hsps[pivot]->StartQ; 
  c1 = i;
  c2 = pivot+1;
  while (c1 < m+i && c2 <= j) 
    if (hsps[c1]->StartQ <= pivotStartQ)
	c1++;
    else
      {
	tmp = hsps[c1];
	hsps[c1] = hsps[c2];
	hsps[c2] = tmp;

	c2++;
      }

  return(pivot);
}

void quickSort(hsp** hsps, int i, int j)
{
  int pivot;

  if (i == j+1)
    /* Nothing */;
  else
    {
      pivot = Split(hsps,i,j);
      
      quickSort(hsps,i,pivot-1);
      quickSort(hsps,pivot+1,j);
    }
}

void SortHSP(packHSP *allHsp)
{
  int i;

  /* quicksorting all the strands and frames */ 
  for(i=0; i<STRANDS*FRAMES; i++)
    quickSort(allHsp->hsps[i],0,allHsp->nHsps[i]-1);
}

