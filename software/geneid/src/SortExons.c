/*************************************************************************
*                                                                        *
*   Module: SortExons                                                    *
*                                                                        *
*   Sort by acceptor position all predicted exons.                       *
*                                                                        *
*   This file is part of the geneid Distribution                         *
*                                                                        *
*     Copyright (C) 2000 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           * 
*                                                                        *
*  This program is free software; you can redistribute it and/or modify  *
*  it under the terms of the GNU General Public License as published by  *
*  the Free Software Foundation; either version 2 of the License, or     *
*  (at your option) any later version.                                   *
*                                                                        *
*  This program is distributed in the hope that it will be useful,       *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*  GNU General Public License for more details.                          *
*                                                                        *
*  You should have received a copy of the GNU General Public License     *
*  along with this program; if not, write to the Free Software           * 
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             *
*************************************************************************/

/*  $Id: SortExons.c,v 1.3 2001-02-08 09:58:18 eblanco Exp $  */

#include "geneid.h"

extern long NUMEXONS; 
extern int EVD;
extern int FWD,RVS;
extern int scanORF;


struct exonitem 
{
  exonGFF *Exon;
  struct exonitem *nexitem;
} ExonItem;

void FreeItems(struct exonitem *q)
{
  if (q == NULL)
    ;
  else
    {
      FreeItems(q->nexitem);
      free(q);
    }
}

void UpdateList(struct exonitem **p, exonGFF* InputExon)
{
  /* Insert at the end if this list */
  while (*p!=NULL)
    p=&((*p)->nexitem); 
  
  /* New item */
  if ((*p= (struct exonitem *) malloc(sizeof(struct exonitem))) == NULL)
    printError("Not enough space to hold one new exonitem");
  
  (*p)->Exon=InputExon;
  (*p)->nexitem=NULL;
}

void SortExons(packExons* allExons, 
	       packExons* allExons_r, 
	       packEvidence* pv,
	       exonGFF* Exons,	       
	       long l1, long l2)
{ 
  struct exonitem **ExonList, *q;
  long i;
  long acceptor;
  long n;
  int offset;
  long l;
  long left;
  long right;
  char mess[MAXSTRING]; 

  /* Boundaries of sorting-array */
  left = l1 + COFFSET - LENGTHCODON;
  right = l2 + COFFSET;
  l =  right - left + 1;

  /* Allocate memory for array ExonList */
  if ((ExonList = 
       (struct exonitem **) calloc(l, sizeof(struct exonitem *))) == NULL)
    printError("Not enough space to hold exonitems list");

  for (i=0;i<l;i++)
    ExonList[i]=NULL;

  /* Adding forward exons ... */
  for (i=0; i<allExons->nInitialExons; i++) 
    {
      acceptor=(allExons->InitialExons+i)->Acceptor->Position - left;
      (allExons->InitialExons+i)->Strand = '+';
      CorrectExon(allExons->InitialExons+i);
      offset = (allExons->InitialExons+i)->offset1;
      UpdateList(&(ExonList[acceptor+offset]), allExons->InitialExons+i); 
    }

  for (i=0;i<allExons->nInternalExons;i++) 
    {
      acceptor=(allExons->InternalExons+i)->Acceptor->Position - left;
      (allExons->InternalExons+i)->Strand = '+';
      CorrectExon(allExons->InternalExons+i);
      offset = (allExons->InternalExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons->InternalExons+i); 
    }

  for (i=0;i<allExons->nTerminalExons;i++) 
    {
      acceptor=(allExons->TerminalExons+i)->Acceptor->Position - left;
      (allExons->TerminalExons+i)->Strand = '+';
      CorrectExon(allExons->TerminalExons+i);
      offset = (allExons->TerminalExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons->TerminalExons+i); 
    }
  
  for (i=0;i<allExons->nSingles;i++) 
    {
      acceptor=(allExons->Singles+i)->Acceptor->Position - left;
      (allExons->Singles+i)->Strand = '+';
      CorrectExon(allExons->Singles+i);
      offset = (allExons->Singles+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons->Singles+i); 
    }
  
  if (scanORF)
    for (i=0;i<allExons->nORFs;i++) 
      {
	acceptor=(allExons->ORFs+i)->Acceptor->Position - left;
	(allExons->ORFs+i)->Strand = '+';
	CorrectORF(allExons->ORFs+i);
	offset = (allExons->ORFs+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->ORFs+i); 
      }

  /* Adding reverse exons ... */
  for (i=0;i<allExons_r->nInitialExons;i++) 
    {
      acceptor=(allExons_r->InitialExons+i)->Acceptor->Position - left;
      (allExons_r->InitialExons+i)->Strand = '-';
      CorrectExon(allExons_r->InitialExons+i);
      offset = (allExons_r->InitialExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons_r->InitialExons+i); 
    }

  for (i=0;i<allExons_r->nInternalExons;i++) 
    {
      acceptor=(allExons_r->InternalExons+i)->Acceptor->Position - left;
      (allExons_r->InternalExons+i)->Strand = '-';
      CorrectExon(allExons_r->InternalExons+i);
      offset = (allExons_r->InternalExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons_r->InternalExons+i); 
    }

  for (i=0;i<allExons_r->nTerminalExons;i++) 
    {   
      acceptor=(allExons_r->TerminalExons+i)->Acceptor->Position - left;
      (allExons_r->TerminalExons+i)->Strand = '-';
      CorrectExon(allExons_r->TerminalExons+i);
      offset = (allExons_r->TerminalExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons_r->TerminalExons+i); 
    }

  for (i=0;i<allExons_r->nSingles;i++) 
    {
      acceptor=(allExons_r->Singles+i)->Acceptor->Position - left;
      (allExons_r->Singles+i)->Strand = '-';
      CorrectExon(allExons_r->Singles+i);
      offset = (allExons_r->Singles+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons_r->Singles+i); 
    }

  if (scanORF)
    for (i=0;i<allExons_r->nORFs;i++) 
      {
	acceptor=(allExons_r->ORFs+i)->Acceptor->Position - left;
	(allExons_r->ORFs+i)->Strand = '-';
	CorrectORF(allExons_r->ORFs+i);
	offset = (allExons_r->ORFs+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->ORFs+i); 
      }

  /* Adding evidence exons */
  if (EVD)
    for (i = pv->i1vExons; i < pv->i2vExons ; i++) 
      {  
	acceptor=(pv->vExons + i)->Acceptor->Position - left;

	/* Requirement 1: range of values */
	if (acceptor < 0 || acceptor > (l2 + LENGTHCODON))
	  {
	    /* Skip wrong exon (A) */
	    sprintf(mess,"Skipped evidence (range): %ld %ld %c %hd",
		    (pv->vExons + i)->Acceptor->Position,
		    (pv->vExons + i)->Donor->Position,
		    (pv->vExons + i)->Strand,
		    (pv->vExons + i)->Frame);
	    printMess(mess);
	  }
	else
	  {
	    /* Requirement 2: forced strand prediction */
	    if ( (!FWD && (pv->vExons + i)->Strand == '+') ||
		 (!RVS && (pv->vExons + i)->Strand == '-'))
	      {
		/* Skip wrong exon (B) */
		sprintf(mess,"Skipped evidence (strand): %ld %ld %c %hd",
			(pv->vExons + i)->Acceptor->Position,
			(pv->vExons + i)->Donor->Position,
			(pv->vExons + i)->Strand,
			(pv->vExons + i)->Frame);
		printMess(mess);
	      }
	    else
	      UpdateList(&(ExonList[acceptor]), pv->vExons+i);
	  }
      }

  /* Scanning all the table searching the exons(sorted by acceptor) */
  for (i=0, n=0; i<l; i++) 
    {
      q=ExonList[i];
      while (q != NULL)
	{
	  Exons[n].Acceptor = q->Exon->Acceptor;
	  Exons[n].Donor = q->Exon->Donor;
	  strcpy(Exons[n].Type,q->Exon->Type);
	  Exons[n].Frame = q->Exon->Frame;
	  Exons[n].Strand = q->Exon->Strand;
	  Exons[n].Score = q->Exon->Score;
	  Exons[n].PartialScore = q->Exon->PartialScore;
	  Exons[n].SRScore = q->Exon->SRScore;
	  Exons[n].GeneScore = 0; 
	  Exons[n].Remainder = q->Exon->Remainder;
	  Exons[n].Group = q->Exon->Group;
	  Exons[n].evidence = q->Exon->evidence;
	  Exons[n].offset1 = q->Exon->offset1;
	  Exons[n].offset2 = q->Exon->offset2; 
	  Exons[n].selected = 0;

	  n++;
	  
	  if (n == RSORTE*NUMEXONS)
	    printError("Too many predicted exons: Change REXONS parameter");
	  
	  q=q->nexitem;
	}
      /* Free all chained items in q */
      FreeItems(ExonList[i]);       
    }

  /* Free empty array */
  free(ExonList);
}

   
