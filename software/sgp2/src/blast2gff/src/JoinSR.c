/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: JoinHSP.c                                          *
*                                                              *
*   Join consecutive SR with the same score                    *
***************************************************************/

#include "blast2gff.h"

void Merge(sr_t* A, sr_t* B)
{
  /* A will remain and B will be free */
  A->End = B-> End;
  A->next = B->next;

  free(B);
}

void Join1SR(packSR* allSr, int i)
{
  sr_t* A;
  sr_t* B; 
  int stop;

  A = allSr->sr[i];
  
  if (A != NULL)
    {
      stop = (A->next == NULL);
      while(!stop)
	{
	  B = A->next;
	  if (((A->End + 1) == B->Start) && (A->Score == B->Score))
	    {
	      /*  printf("Join %ld %ld -- %ld %ld\n", */
	      /*  		     A->Start,A->End, */
	      /*  		     B->Start,B->End); */
	      Merge(A,B);
	      allSr->nr[i]--;
	    }
	  else
	    {
	      /* Next loop */
	      A = A->next;
	    }
	  stop = (A->next == NULL);
	}
    }
}

void JoinSR(packSR* allSr)
{
  int i;

  /* join consecutive SR with the same score */
  for(i=0; i<STRANDS*FRAMES; i++)
    Join1SR(allSr,i);  
}
