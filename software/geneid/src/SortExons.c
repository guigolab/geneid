/*************************************************************************
*                                                                        *
*   Module: SortExons                                                    *
*                                                                        *
*   Sort by left signal position all of predicted exons                  *
*                                                                        *
*   This file is part of the geneid 1.2 distribution                     *
*                                                                        *
*     Copyright (C) 2003 - Enrique BLANCO GARCIA                         *
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

/*  $Id: SortExons.c,v 1.6 2003-11-05 15:09:25 eblanco Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic exons (multiplied by FSORT) */
extern long NUMEXONS;

extern int EVD;
extern int FWD,RVS;
extern int scanORF;
extern int SGE;

/* Artificial initial gene feature: force complete gene prediction */
void InsertFirstGhostExon(exonGFF* Exons)
{
  exonGFF* e;
  exonGFF* re;

  /* 1. Forward Strand */
  /* Create the exon structure */
  if ((e = (exonGFF*) malloc(sizeof(exonGFF))) == NULL)
        printError("Not enough memory: first ghost exon(sGHOST)");

  /* Create the sites structure */
  if ((e->Acceptor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor site for first ghost exon(sGHOST)");

  /* Create the sites structure */
  if ((e->Donor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor site for first ghost exon(sGHOST)");
  
  e->Acceptor->Position = 0;
  e->Donor->Position = 0;
  
  /* Save the extracted exon */
  Exons[0].Acceptor = e->Acceptor;
  Exons[0].Donor = e->Donor;
  strcpy(Exons[0].Type,sGHOST);
  Exons[0].Frame = 0;
  Exons[0].Strand = '+';
  Exons[0].Score = MAXSCORE;
  Exons[0].PartialScore = 0;
  Exons[0].HSPScore = 0;
  Exons[0].GeneScore = 0; 
  Exons[0].Remainder = 0;
  strcpy(Exons[0].Group,NOGROUP);
  Exons[0].evidence = 0;
  Exons[0].offset1 = 0;
  Exons[0].offset2 = 0; 
  Exons[0].lValue = 0;
  Exons[0].rValue = 0; 
  Exons[0].selected = 0;

  /* 2. Reverse Strand */
  /* Create the exon structure */
  if ((re = (exonGFF*) malloc(sizeof(exonGFF))) == NULL)
        printError("Not enough memory: first ghost exon(sGHOST)");

  /* Create the sites structure */
  if ((re->Acceptor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor site for first ghost exon(sGHOST)");

  /* Create the sites structure */
  if ((re->Donor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor site for first ghost exon(sGHOST)");
  
  re->Acceptor->Position = 0;
  re->Donor->Position = 0;
  
  /* Save the extracted exon */
  Exons[1].Acceptor = re->Acceptor;
  Exons[1].Donor = re->Donor;
  strcpy(Exons[1].Type,sGHOST);
  Exons[1].Frame = 0;
  Exons[1].Strand = '-';
  Exons[1].Score = MAXSCORE;
  Exons[1].PartialScore = 0;
  Exons[1].HSPScore = 0;
  Exons[1].GeneScore = 0; 
  Exons[1].Remainder = 0;
  strcpy(Exons[1].Group,NOGROUP);
  Exons[1].evidence = 0;
  Exons[1].offset1 = 0;
  Exons[1].offset2 = 0; 
  Exons[1].lValue = 0;
  Exons[1].rValue = 0; 
  Exons[1].selected = 0;
}

/* Artificial terminal gene feature: force complete gene prediction */
void InsertLastGhostExon(exonGFF* Exons, long n, long L)
{
  exonGFF* e;
  exonGFF* re;

  /* 1. Forward Strand */
  /* Create the exon structure */
  if ((e = (exonGFF*) malloc(sizeof(exonGFF))) == NULL)
        printError("Not enough memory: last ghost exon(sGHOST)");

  /* Create the sites structure */
  if ((e->Acceptor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor site for last ghost exon(sGHOST)");

  /* Create the sites structure */
  if ((e->Donor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor site for last ghost exon(sGHOST)");

  e->Acceptor->Position = L;
  e->Donor->Position = L;

  /* Save the extracted exon */
  Exons[n].Acceptor = e->Acceptor;
  Exons[n].Donor = e->Donor;
  strcpy(Exons[n].Type,sGHOST);
  Exons[n].Frame = 0;
  Exons[n].Strand = '+';
  Exons[n].Score = MAXSCORE;
  Exons[n].PartialScore = 0;
  Exons[n].HSPScore = 0;
  Exons[n].GeneScore = 0; 
  Exons[n].Remainder = 0;
  strcpy(Exons[n].Group,NOGROUP);
  Exons[n].evidence = 0;
  Exons[n].offset1 = 0;
  Exons[n].offset2 = 0; 
  Exons[n].lValue = 0;
  Exons[n].rValue = 0; 
  Exons[n].selected = 0;

  /* 1. Reverse Strand */
  /* Create the exon structure */
  if ((re = (exonGFF*) malloc(sizeof(exonGFF))) == NULL)
        printError("Not enough memory: last ghost exon(sGHOST)");

  /* Create the sites structure */
  if ((re->Acceptor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor site for last ghost exon(sGHOST)");

  /* Create the sites structure */
  if ((re->Donor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor site for last ghost exon(sGHOST)");

  re->Acceptor->Position = L;
  re->Donor->Position = L;

  /* Save the extracted exon */
  Exons[n+1].Acceptor = re->Acceptor;
  Exons[n+1].Donor = re->Donor;
  strcpy(Exons[n+1].Type,sGHOST);
  Exons[n+1].Frame = 0;
  Exons[n+1].Strand = '-';
  Exons[n+1].Score = MAXSCORE;
  Exons[n+1].PartialScore = 0;
  Exons[n+1].HSPScore = 0;
  Exons[n+1].GeneScore = 0; 
  Exons[n+1].Remainder = 0;
  strcpy(Exons[n+1].Group,NOGROUP);
  Exons[n+1].evidence = 0;
  Exons[n+1].offset1 = 0;
  Exons[n+1].offset2 = 0; 
  Exons[n+1].lValue = 0;
  Exons[n+1].rValue = 0; 
  Exons[n+1].selected = 0;
}

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
               packExternalInformation* external,
			   packEvidence* pv,
			   exonGFF* Exons,         
               long l1, long l2,
			   long LengthSequence)
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
  if (EVD && pv != NULL)
    for (i = external->i1vExons; i < external->i2vExons ; i++) 
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


  /* SGE- STARTING fragment: Insert first ghost exon sGHOST */
  if (SGE && (l1 == 0))
	{
	  n = 2;
	  InsertFirstGhostExon(Exons);
	}
  else
	n = 0;

  /* 3. Traversing the table extracting the exons sorted by left position */
  HowMany = FSORT*NUMEXONS;
  /* for (i=0, n=0; i<l; i++) */
  for (i=0; i<l; i++)
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
		  Exons[n].HSPScore = q->Exon->HSPScore;
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

  /* SGE- FINISHING split: Insert last ghost exon "$" */
  if (SGE && (l2 == LengthSequence - 1))
	{
	  InsertLastGhostExon(Exons, n, LengthSequence);
	  n = n + 2;
	  
	  if (n >= HowMany)
		printError("Too many predicted exons: increase FSORT parameter");
	}
  
  /* Free empty array */
  free(ExonList);
}

   
