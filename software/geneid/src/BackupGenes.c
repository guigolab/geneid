/*************************************************************************
*                                                                        *
*   Module: BackupGenes                                                  *
*                                                                        *
*   To save best partial genes between 2 contigous fragments             *
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

/*  $Id: BackupGenes.c,v 1.4 2003-11-05 13:23:01 eblanco Exp $  */

#include "geneid.h"

/* Maximum allowed number of sites and exons to save */
extern long MAXBACKUPSITES, MAXBACKUPEXONS;

/* Increase counters modulus a long number */
long IncrMod(long x, long Modulus)
{
  long z;
  
  z = x+1;
  if (z == Modulus)
    z = 0;
  
  return(z);
}

/* Saving exon information (features) into the dumpster */
exonGFF* backupExon(exonGFF* E, exonGFF* Prev, packDump* d)
{
  /* back-up acceptor */
  d->dumpSites[d->ndumpSites].Position = E->Acceptor->Position;
  d->dumpSites[d->ndumpSites].Score = E->Acceptor->Score;
  d->dumpExons[d->ndumpExons].Acceptor = &(d->dumpSites[d->ndumpSites]); 
  d->ndumpSites = IncrMod(d->ndumpSites, MAXBACKUPSITES);
  
  /* back-up donor */
  d->dumpSites[d->ndumpSites].Position = E->Donor->Position;
  d->dumpSites[d->ndumpSites].Score = E->Donor->Score;
  d->dumpExons[d->ndumpExons].Donor = &(d->dumpSites[d->ndumpSites]); 
  d->ndumpSites = IncrMod(d->ndumpSites, MAXBACKUPSITES);
  
  /* back-up exon properties */
  strcpy(d->dumpExons[d->ndumpExons].Type, E->Type);
  d->dumpExons[d->ndumpExons].Frame  = E->Frame;
  d->dumpExons[d->ndumpExons].Strand = E->Strand; 
  d->dumpExons[d->ndumpExons].Score  = E->Score;
  d->dumpExons[d->ndumpExons].PartialScore = E->PartialScore;
  d->dumpExons[d->ndumpExons].HSPScore = E->HSPScore;
  d->dumpExons[d->ndumpExons].GeneScore  = E->GeneScore;
  d->dumpExons[d->ndumpExons].Remainder = E->Remainder; 
  strcpy(d->dumpExons[d->ndumpExons].Group,E->Group);
  d->dumpExons[d->ndumpExons].offset1 = E->offset1;
  d->dumpExons[d->ndumpExons].offset2 = E->offset2;
  d->dumpExons[d->ndumpExons].lValue = E->lValue;
  d->dumpExons[d->ndumpExons].rValue = E->rValue;
  d->dumpExons[d->ndumpExons].evidence = E->evidence;
  d->dumpExons[d->ndumpExons].PreviousExon = Prev;
  
  /* Returns the new exon recently created */
  return(&(d->dumpExons[d->ndumpExons]));
}

/* Saving all about a gene: exons, sites, properties */
exonGFF* backupGene(exonGFF* E, packDump* d)
{
  exonGFF* PrevExon;
  exonGFF* ResExon;
  
  /* Ghost exon doesn't need backup */
  if (E->Strand == '*')
    ResExon = E; 
  else
	{
	  /* Ckeckpoint to discover if exon is already in the dumpster */
	  ResExon  = (exonGFF*) getExonDumpHash(E, d->h);
	  
	  /* New exon: save it and insert into the hash table */
	  if (ResExon == NULL)
		{
		  PrevExon = backupGene(E->PreviousExon,d);
		  ResExon = backupExon(E,PrevExon,d);
		  d->ndumpExons = IncrMod(d->ndumpExons, MAXBACKUPEXONS);
		  /* adding this exon at hash table */
		  setExonDumpHash(ResExon, d->h);       
		}
	  /* if this exon exists, finish backup gene */
	}
  return(ResExon);
}

/* It saves the information about partial genes (packGenes) */
void BackupGenes(packGenes* pg, int nclass, packDump* d)
{
  int i,j;
  
  /* 1. back-up best partial genes */
  for(i=0; i<nclass; i++)
    for(j=0; j<FRAMES; j++)
      pg->Ga[i][j] = backupGene(pg->Ga[i][j], d);
  
  /* 2. back-up optimal(partial gene) */
  pg->GOptim = backupGene(pg->GOptim, d);
}

/* It saves information about d-exons: exons needed the next iteration */
void BackupArrayD(packGenes* pg, long accSearch,
                  gparam* gp, packDump* dumpster)
{
  int i;
  long j;
  long jUpdate,jMaxdist;
  long MinDist;
  long MaxDist;
  long nBackups=0;
  short remainder;
  char mess[MAXSTRING];
  
  /* Traversing sort-by-donor array to save some genes (assembling rules) */
  /* These exons have to be beyond the point accSearch (preserve ordering) */
  for(i=0; i < gp->nclass ; i++)
    {
      /* Get data from this assembling rule */
      j = pg->je[i];
      MaxDist = gp->Md[i];
      MinDist = gp->md[i];
      
      /* 1. Update best genes using remaining exons before accSearch - md */
      while(j < pg->km[i] &&
			((pg->d[i][j]->Donor->Position + pg->d[i][j]->offset2)
			 <=
			 accSearch - MinDist))
		{
		  if (pg->d[i][j]->Strand == cFORWARD)
			remainder = pg->d[i][j]->Remainder;
		  else
			remainder = pg->d[i][j]->Frame;
		  
		  if (pg->d[i][j]->GeneScore > pg->Ga[i][remainder]->GeneScore)
			pg->Ga[i][remainder] = pg->d[i][j];
		  j++;
		}
      jUpdate = j;
      
      /* 2. Copy exons needed to recompute maximum distance requirement */
      if (MaxDist == INFI)      
		jMaxdist = jUpdate;
      else
		{
		  /* Loop back until reaching first exon out of range: acc-dMax */
		  j = jUpdate;
		  while (j>=0 && j < pg->km[i]  && 
				 ((pg->d[i][j]->Donor->Position 
				   + pg->d[i][j]->offset2)
				  >= 
				  accSearch - MaxDist))
			j--;
		  jMaxdist = (j < pg->km[i])? j+1 : j;
		}
      
      /* 3. To save the set of best genes finished in any selected d-exon */
      for(j = jMaxdist; j < pg->km[i]; j++)
		{
		  pg->d[i][j-jMaxdist] =  backupGene(pg->d[i][j], dumpster);
		  nBackups++;
		}
      pg->km[i] = pg->km[i] - jMaxdist;
      pg->je[i] = jUpdate - jMaxdist;
    }
  
  sprintf(mess,"%ld d-genes saved(%ld real exons)",
		  nBackups, dumpster->h->total);
  printMess(mess);
}

/* Reset counters and pointers for the next input sequence */
void cleanGenes(packGenes* pg, int nclass, packDump* dumpster)
{
  int aux, aux2;
  
  for(aux=0; aux<nclass; aux++)
	{
      /* Reset sort-by-donor functions */
      pg->je[aux] = 0;
      pg->km[aux] = 0;
      
      /* Reset Ga-exons: every Ga looks at Ghost exon */
      for(aux2=0; aux2 < FRAMES; aux2++)
		pg->Ga[aux][aux2] = pg->Ghost;
	}
  
  /* Reset Optimal Gene */
  pg->GOptim = pg->Ghost;
  
  /* Reset counters for dumpster arrays if were used before */
  if (MAXBACKUPSITES && MAXBACKUPEXONS)
	{
	  dumpster->ndumpExons = 0;
	  dumpster->ndumpSites = 0;
	}
}
