/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: Output                                             *
*                                                              *
*   Management of options for printing results                 *
***************************************************************/   

#include "blast2gff.h"

extern int VRB, GFFOUT, PSR;

extern account *m;    

/* Printing messages (information) */
void printMess(char* s)
{
  if (VRB)
     fprintf (stderr, "> %s\n",s);
}

/* Printing partial and final results (information) */
void printRes(char* s)
{
  if (VRB)
     fprintf (stderr, "\t%s\n",s);
}

/* Printing error messages */
void printError(char *s)
{
     fprintf(stderr,"Error: %s\n",s);
	  exit(1);
}

/* Get time of execution using account information */
void OutputTime()
{
  time_t tEnd;
  int t;
  float caux;
  char mess[MAXSTRING];

  (void) time(&tEnd);
  /* caux = (float)clock() / (float)CLOCKS_PER_SEC; */
  caux = 0;
  t = (int) tEnd - m->tStart;

  if (t < caux)
    t++;

  printRes("_______________________________________________________________\n");
 
  sprintf(mess,"CPU time: \t%.3f secs",caux);
  printRes(mess);

  sprintf(mess,"Total time: \t%d secs(%.2f mins)\n\n",t,(float)t/MINUTE);
  printRes(mess);
}

void PrintHSP(packHSP* allHSP, char Query[LOCUSLENGTH])
{
  int i;
  long j;

  for(i=0; i < FRAMES*STRANDS; i++)
    for(j=0; j < allHSP->nHsps[i]; j++)
      printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\t%s\n",
	     Query,
	     VERSION,
	     HSP,
	     allHSP->hsps[i][j]->StartQ,
	     allHSP->hsps[i][j]->EndQ,
	     allHSP->hsps[i][j]->Score,
	     allHSP->hsps[i][j]->StrandQ,
	     allHSP->hsps[i][j]->FrameQ,
	     allHSP->hsps[i][j]->Subject);
}

void PrintSRStrands(packSR* allSR)
{
  int i;
  sr_t* tmp;

  printf("## blast2gff v 1.0. imim.es\n");

  for(i=0; i< FRAMES*STRANDS; i++)
    {
      printf("## Created %d SR(frame %d)\n",allSR->nr[i],i);

      tmp = allSR->sr[i];
      while(tmp != NULL)
	{
	  printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\n",
		 NONAME,
		 VERSION,
		 SR,
		 tmp->Start,
		 tmp->End,
		 tmp->Score,
		 tmp->Strand,
		 tmp->Frame); 
	  tmp = tmp->next;
	}
    }
}

void PrintSR(packSR* allSR, char Query[LOCUSLENGTH])
{
  int i;

  printf("## blast2gff v 1.0. imim.es\n");
  printf("## Created %ld SR\n",allSR->nSR);

  for(i=0; i < allSR->nSR; i++)
    {
      printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\n",
	     Query,
	     VERSION,
	     SR,
	     allSR->sortSR[i]->Start,
	     allSR->sortSR[i]->End,
	     allSR->sortSR[i]->Score,
	     allSR->sortSR[i]->Strand,
	     allSR->sortSR[i]->Frame); 
    }
}


void Output(packHSP* allHSP, packSR* allSR, char Query[LOCUSLENGTH])
{
  if (GFFOUT)
    PrintHSP(allHSP,Query);

  if (PSR)
      PrintSR(allSR,Query);
}






