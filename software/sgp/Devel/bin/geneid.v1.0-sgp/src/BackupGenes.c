/*************************************************************************
*                                                                        *
*   Module: BackupGenes                                                  *
*                                                                        *
*   Save genes that GenAmic will need at next split.                     *
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

/*  $Id: BackupGenes.c,v 1.2 2000-07-26 08:25:19 eblanco Exp $  */

#include "geneid.h"

extern long MAXBACKUPSITES, MAXBACKUPEXONS;

long IncrMod(long x, long Modulus)
{
  long z;
  
  z = x+1;
  if (z == Modulus)
    z = 0;
  return(z);
}

/* It saves the information about a exon */
exonGFF* backupExon(exonGFF* E, exonGFF* Prev, packDump* d)
{
  /* backup acceptor */
  d->dumpSites[d->ndumpSites].Position = E->Acceptor->Position;
  d->dumpSites[d->ndumpSites].Score = E->Acceptor->Score;
  d->dumpExons[d->ndumpExons].Acceptor = &(d->dumpSites[d->ndumpSites]); 
  d->ndumpSites = IncrMod(d->ndumpSites, MAXBACKUPSITES);
  
  /* backup donor */
  d->dumpSites[d->ndumpSites].Position = E->Donor->Position;
  d->dumpSites[d->ndumpSites].Score = E->Donor->Score;
  d->dumpExons[d->ndumpExons].Donor = &(d->dumpSites[d->ndumpSites]); 
  d->ndumpSites = IncrMod(d->ndumpSites, MAXBACKUPSITES);
  
  /* back-up exon */
  strcpy(d->dumpExons[d->ndumpExons].Type, E->Type);
  d->dumpExons[d->ndumpExons].Frame  = E->Frame;
  d->dumpExons[d->ndumpExons].Strand = E->Strand; 
  d->dumpExons[d->ndumpExons].Score  = E->Score;
  d->dumpExons[d->ndumpExons].PartialScore = E->PartialScore;
  d->dumpExons[d->ndumpExons].SRScore = E->SRScore;
  d->dumpExons[d->ndumpExons].GeneScore  = E->GeneScore;
  d->dumpExons[d->ndumpExons].Remainder = E->Remainder; 
  d->dumpExons[d->ndumpExons].Group = E->Group;
  d->dumpExons[d->ndumpExons].offset1 = E->offset1;
  d->dumpExons[d->ndumpExons].offset2 = E->offset2;
  d->dumpExons[d->ndumpExons].evidence = E->evidence;
  d->dumpExons[d->ndumpExons].PreviousExon = Prev;
  
  return(&(d->dumpExons[d->ndumpExons]));
}

/* It saves the information about a gene, saving all its exons */
exonGFF* backupGene(exonGFF* E, packDump* d)
{
  exonGFF* PrevExon;
  exonGFF* ResExon;
  
  /* Ghost exon doesn't need backup */
  if (E->Strand == '*')
    ResExon = E; 
  else
     {
       /* Is this exon backuped before? */
       ResExon  = (exonGFF*) getExonDumpHash(E, d->h);
       
       /* New exon: backup it and add it at hash table */
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

/* It saves the information of partial genes */
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

/* It saves the information of d-genes */
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

  for(i=0; i < gp->nclass ; i++)
    {
      j = pg->je[i];
      MaxDist = gp->Md[i];
      MinDist = gp->md[i];
      
      /* First: update d-array until the boundary-exon acceptor */
      while(j < pg->km[i] &&
	    ((pg->d[i][j]->Donor->Position + pg->d[i][j]->offset2)
	     <=
	     accSearch - MinDist))
	{
	  remainder = pg->d[i][j]->Remainder;
	  if (pg->d[i][j]->GeneScore > pg->Ga[i][remainder]->GeneScore)
	    pg->Ga[i][remainder] = pg->d[i][j];
	  j++;
	}
      jUpdate = j;
      
      /* Second: How many d-exons must we backup? */
      if (MaxDist == INFI)      
	jMaxdist = jUpdate;
      else
	{
	  /* Loop back until first exon out of acc-dMax */
	  j = jUpdate;
	  while (j>=0 && j < pg->km[i]  && 
		 ((pg->d[i][j]->Donor->Position 
		   + pg->d[i][j]->offset2)
		  >= 
		  accSearch - MaxDist))
            j--;
	  jMaxdist = (j < pg->km[i])? j+1 : j;
	}
      
      /* Third: save the current state of d-array */
      for(j = jMaxdist; j < pg->km[i]; j++)
	{
	  pg->d[i][j-jMaxdist] =  backupGene(pg->d[i][j], dumpster);
	  nBackups++;
	}
      pg->km[i] = pg->km[i] - jMaxdist;
      pg->je[i] = jUpdate - jMaxdist;
    }
 
  sprintf(mess,"%ld d-genes saved(%ld real exons)\n",
	  nBackups, dumpster->h->total);
  printMess(mess);

}

/* Clean old backup genes because next sequence begins */
void cleanGenes(packGenes* pg, int nclass, packDump* dumpster)
{
   int aux, aux2;

  for(aux=0; aux<nclass; aux++)
  {
      /* Clean all sort-by-donor array */
      pg->je[aux] = 0;
      pg->km[aux] = 0;
      
      /* Reinitialize Ga-exons: every Ga looks at Ghost exon */
      for(aux2=0; aux2 < FRAMES; aux2++)
	pg->Ga[aux][aux2] = pg->Ghost;
  }

  /* Clean Optimal Gen */
  pg->GOptim = pg->Ghost;

  if (MAXBACKUPSITES && MAXBACKUPEXONS)
  {
    /* Reset dumpster */
    dumpster->ndumpExons = 0;
    dumpster->ndumpSites = 0;
  }
}




