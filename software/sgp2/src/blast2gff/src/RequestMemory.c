/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: RequestMemory                                      *
*                                                              *
*    Request memory to System to store Blast2GFF data structs  *
***************************************************************/

#include "blast2gff.h"

packHSP* RequestMemoryHSP()
{
  packHSP *allHsp;
  int i; 

  if ((allHsp =
       (packHSP*) malloc(sizeof(packHSP)))  == NULL)
    printError("Not enough space to hold HSP data structure");

  /* HSP */
  if ((allHsp->hsps =
       (hsp ** *) calloc(STRANDS*FRAMES, sizeof(hsp**))) == NULL)
    printError("Not enough space to hold HSP 6-array");  


  for(i=0;i<STRANDS*FRAMES;i++)
    if ((allHsp->hsps[i] =
	 (hsp**) calloc(NUMHSPS, sizeof(hsp*)))  == NULL)
      printError("Not enough space to hold HSPs");

  /* Counters */
  if ((allHsp->nHsps =
       (long*) calloc(STRANDS*FRAMES, sizeof(int)))  == NULL)
    printError("Not enough space to hold HSP counters");  

  /* Hack of Alpha */
  for(i=0;i<STRANDS*FRAMES;i++)
    allHsp->nHsps[i] = 0;

  allHsp->nTotalHsps = 0;

  return(allHsp);
}

packSR* RequestMemorySR()
{
  packSR *allSr;
  int i;

  if ((allSr =
       (packSR*) malloc(sizeof(packSR)))  == NULL)
    printError("Not enough space to hold SR data structure");

  for(i=0;i<STRANDS*FRAMES;i++)
    {
      allSr->nr[i]  = 0;
      allSr->nr0[i] = NULL;
      allSr->sr[i]  = NULL;
    }
  
  allSr->nSR = 0;
  
  return(allSr);
}
