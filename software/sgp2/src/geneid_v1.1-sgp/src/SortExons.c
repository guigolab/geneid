/*************************************************************************
*                                                                        *
*   Module: SortExons                                                    *
*                                                                        *
*   Sort by left signal position all of predicted exons                  *
*                                                                        *
*   This file is part of the geneid 1.1 distribution                     *
*                                                                        *
*     Copyright (C) 2001 - Enrique BLANCO GARCIA                         *
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

/*  $Id: SortExons.c,v 1.1 2003-09-10 14:53:34 gparra Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic exons (multiplied by FSORT) */
extern long NUMEXONS;

extern int EVD;
extern int FWD,RVS;
extern int scanORF;

/* Struct for a node (list): pointer to exon and to next node */
struct exonitem 
{
  exonGFF* Exon;
  struct exonitem* nexitem;
} ExonItem;

/* Set free a list of nodes */
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

/* Insert an exon in the selected list according to its left position */
void UpdateList(struct exonitem** p, exonGFF* InputExon)
{
  /* Insert the new node at the end of the list */
  while (*p!=NULL)
    p=&((*p)->nexitem); 
  
  /* Allocating new node for this exon */
  if ((*p= (struct exonitem *) malloc(sizeof(struct exonitem))) == NULL)
    printError("Not enough memory: new exonitem (sorting)");
  
  /* Updating information for the node */
  (*p)->Exon=InputExon;
  (*p)->nexitem=NULL;
}

/* Sort all of predicted exons by placing them in an array of lists */
/* corresponding every list to a beginning position for predicted exons */
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
  long room;
  char mess[MAXSTRING]; 
  long HowMany;

  /* 0. Creating the array for sorting: 1 - Length of fragment */
  left = l1 + COFFSET - LENGTHCODON;
  right = l2 + COFFSET;
  l =  right - left + 1;

  /* Allocating memory for array ExonList */
  if ((ExonList = 
       (struct exonitem **) calloc(l, sizeof(struct exonitem *))) == NULL)
    printError("Not enough memory: ExonList array (sorting)");

  /* Reset the positions, pointing to NULL */
  for (i=0;i<l;i++)
    ExonList[i]=NULL;

  /* 1. Insert exons in the proper list according to its left position */
  /* Adding predicted exons in the forward sense */
  for (i=0; i<allExons->nInitialExons; i++) 
    {
      /* Correction between real and relative to fragment position */
      acceptor=(allExons->InitialExons+i)->Acceptor->Position - left;
      (allExons->InitialExons+i)->Strand = '+';
      /* Fixing positions according to the COFFSET in arrays */
      CorrectExon(allExons->InitialExons+i);
      offset = (allExons->InitialExons+i)->offset1;
      /* Insert exon in the proper list depending on the left position */
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

  /* Adding predicted exons in the reverse sense */
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

  /* 2. Merge evidence exons (annotations) with predictions */
  if (EVD)
    for (i = pv->i1vExons; i < pv->i2vExons ; i++) 
      {  
		acceptor=(pv->vExons + i)->Acceptor->Position - left;
		offset = (pv->vExons + i)->offset1;
		room = ((acceptor + offset)<0)? 0 : (acceptor + offset);

		/* Requirement 1: range of values */
		if (acceptor < 0)
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
			  UpdateList(&(ExonList[room]), pv->vExons+i);
		  }
      }

  /* 3. Traversing the table extracting the exons sorted by left position */
  HowMany = FSORT*NUMEXONS;
  for (i=0, n=0; i<l; i++)
    {
      /* Extracting exons from list q */
      q=ExonList[i];
      while (q != NULL)
		{
		  /* Save the extracted exon */
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
		  strcpy(Exons[n].Group,q->Exon->Group);
		  Exons[n].evidence = q->Exon->evidence;
		  Exons[n].offset1 = q->Exon->offset1;
		  Exons[n].offset2 = q->Exon->offset2; 
		  Exons[n].lValue = q->Exon->lValue;
		  Exons[n].rValue = q->Exon->rValue; 
		  Exons[n].selected = 0;

		  n++;
	  
		  if (n >= HowMany)
			printError("Too many predicted exons: increase FSORT parameter");
	  
		  q=q->nexitem;
		}

	  /* Free chained items in the processed list q */
	  FreeItems(ExonList[i]);       
	}

  /* Free empty array */
  free(ExonList);
}

   
