/*************************************************************************
*                                                                        *
*   Module: BackupGenes                                                  *
*                                                                        *
*   To save best partial genes between 2 contigous fragments             *
*                                                                        *
*   This file is part of the geneid 1.4 distribution                     *
*                                                                        *
*     Copyright (C) 2006 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
*                          Tyler   ALIOTO                                *
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

/*  $Id: BackupGenes.c,v 1.11 2011-01-13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum allowed number of sites and exons to save */
extern long MAXBACKUPSITES, MAXBACKUPEXONS;
/* The number of compatible splice classes */
extern short SPLICECLASSES;

/* Return a stable-address pointer to exon backup slot i, allocating its chunk
   on first touch. Chunks are never moved (only the array of chunk pointers may
   realloc), so a pointer handed out here stays valid for the whole run. */
static exonGFF* dumpExonAt(packDump* d, long i)
{
  long c = i / DUMPCHUNK;

  if (c >= d->dumpExonsChunks)
    {
      long ncap = d->dumpExonsChunks ? d->dumpExonsChunks : 1;
      exonGFF** tmp;
      long k;
      while (ncap <= c)
	ncap *= 2;
      if ((tmp = (exonGFF**) realloc(d->dumpExons, ncap * sizeof(exonGFF*))) == NULL)
	printError("Not enough memory: backup exon chunk table");
      for (k = d->dumpExonsChunks; k < ncap; k++)
	tmp[k] = NULL;
      d->dumpExons = tmp;
      d->dumpExonsChunks = ncap;
    }
  if (d->dumpExons[c] == NULL)
    if ((d->dumpExons[c] = (exonGFF*) calloc(DUMPCHUNK, sizeof(exonGFF))) == NULL)
      printError("Not enough memory: backup exon chunk");

  return &(d->dumpExons[c][i % DUMPCHUNK]);
}

/* As dumpExonAt, for the site backup array. */
static site* dumpSiteAt(packDump* d, long i)
{
  long c = i / DUMPCHUNK;

  if (c >= d->dumpSitesChunks)
    {
      long ncap = d->dumpSitesChunks ? d->dumpSitesChunks : 1;
      site** tmp;
      long k;
      while (ncap <= c)
	ncap *= 2;
      if ((tmp = (site**) realloc(d->dumpSites, ncap * sizeof(site*))) == NULL)
	printError("Not enough memory: backup site chunk table");
      for (k = d->dumpSitesChunks; k < ncap; k++)
	tmp[k] = NULL;
      d->dumpSites = tmp;
      d->dumpSitesChunks = ncap;
    }
  if (d->dumpSites[c] == NULL)
    if ((d->dumpSites[c] = (site*) calloc(DUMPCHUNK, sizeof(site))) == NULL)
      printError("Not enough memory: backup site chunk");

  return &(d->dumpSites[c][i % DUMPCHUNK]);
}

/* Saving exon information (features) into the dumpster */
exonGFF* backupExon(exonGFF* E, exonGFF* Prev, packDump* d)
{
  exonGFF* de = dumpExonAt(d, d->ndumpExons);
  site* ds;

  /* back-up acceptor */
  ds = dumpSiteAt(d, d->ndumpSites);
  ds->Position = E->Acceptor->Position;
  ds->Score = E->Acceptor->Score;
  ds->ScoreAccProfile = E->Acceptor->ScoreAccProfile;
  ds->ScorePPT = E->Acceptor->ScorePPT;
  ds->ScoreBP = E->Acceptor->ScoreBP;
  ds->PositionBP = E->Acceptor->PositionBP;
  ds->PositionPPT = E->Acceptor->PositionPPT;
  ds->class = E->Acceptor->class;
  strcpy(ds->subtype,E->Acceptor->subtype);
  strcpy(ds->type,E->Acceptor->type);
  de->Acceptor = ds;
  d->ndumpSites = d->ndumpSites + 1;

  /* back-up donor */
  ds = dumpSiteAt(d, d->ndumpSites);
  ds->Position = E->Donor->Position;
  ds->Score = E->Donor->Score;
  ds->ScoreAccProfile = E->Donor->ScoreAccProfile;
  ds->ScorePPT = E->Donor->ScorePPT;
  ds->ScoreBP = E->Donor->ScoreBP;
  ds->PositionBP = E->Donor->PositionBP;
  ds->PositionPPT = E->Donor->PositionPPT;
  ds->class = E->Donor->class;
  strcpy(ds->subtype,E->Donor->subtype);
  strcpy(ds->type,E->Donor->type);
  de->Donor = ds;
  d->ndumpSites = d->ndumpSites + 1;

  /* back-up exon properties */
  strcpy(de->Type, E->Type);
  de->Frame  = E->Frame;
  de->Strand = E->Strand;
  de->Score  = E->Score;
  de->PartialScore = E->PartialScore;
  de->HSPScore = E->HSPScore;
  de->R = E->R;
  de->GeneScore  = E->GeneScore;
  de->Remainder = E->Remainder;
  strcpy(de->Group,E->Group);
  de->offset1 = E->offset1;
  de->offset2 = E->offset2;
  de->lValue = E->lValue;
  de->rValue = E->rValue;
  de->evidence = E->evidence;
  de->PreviousExon = Prev;

  /* Returns the new exon recently created */
  return(de);
}

/* Saving all about a gene: exons, sites, properties */
exonGFF* backupGene(exonGFF* E, packDump* d)
{
  exonGFF* PrevExon;
  exonGFF* ResExon;

  /* Ghost exon doesn't need backup */
  if ((E->Strand == '*')) /* ||(E->Strand != '+')||(E->Strand != '-')) */
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


		  d->ndumpExons = d->ndumpExons + 1;
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
  int i,j,k;

  /* 1. back-up best partial genes */
  for(i=0; i<nclass; i++)
    for(j=0; j<FRAMES; j++)
      for(k=0; k<SPLICECLASSES; k++){
	pg->Ga[i][j][k] = backupGene(pg->Ga[i][j][k], d);
      }

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
  short donorclass;
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
		  if (pg->d[i][j]->Strand == cFORWARD){
			remainder = pg->d[i][j]->Remainder;
			donorclass = pg->d[i][j]->Donor->class;
		  }else{
			remainder = pg->d[i][j]->Frame;
			donorclass = pg->d[i][j]->Acceptor->class;
		  }
		  
		  if (pg->d[i][j]->GeneScore > pg->Ga[i][remainder][donorclass]->GeneScore)
			pg->Ga[i][remainder][donorclass] = pg->d[i][j];
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
  int aux, aux2, aux3;
/*   for(aux=0; aux<nclass; aux++) */
  for(aux=0; aux<nclass; aux++)
    {
      /* Reset sort-by-donor functions */
      pg->je[aux] = 0;
      pg->km[aux] = 0;
      
      /* Reset Ga-exons: every Ga looks at Ghost exon */
      for(aux2=0; aux2 < FRAMES; aux2++){
	for(aux3=0; aux3 < SPLICECLASSES; aux3++){
	  pg->Ga[aux][aux2][aux3] = pg->Ghost;
	}
      }
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
