/*************************************************************************
*                                                                        *
*   Module: RequestMemory                                                *
*                                                                        *
*   Request memory to System to store geneid data structures             *
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

/*  $Id: RequestMemory.c,v 1.2 2000-07-26 07:59:10 eblanco Exp $  */

#include "geneid.h"

extern long NUMSITES, NUMEXONS, MAXBACKUPSITES, MAXBACKUPEXONS;

void shareGeneModel(gparam** isochores, int n)
{ 
  int i,j,k;
 
  /* Sharing global parameters */
  for(i=1; i<n; i++)
    {
      isochores[i]->D      = isochores[0]->D; 
      isochores[i]->nc     = isochores[0]->nc; 
      isochores[i]->ne     = isochores[0]->ne; 
      isochores[i]->md     = isochores[0]->md;
      isochores[i]->Md     = isochores[0]->Md;
      isochores[i]->nclass = isochores[0]->nclass;

      /* Copy class info */
      for(j=0; j<isochores[i]->nclass; j++)
	{
	  for(k=0; k< isochores[i]->nc[j]; k++)
	    isochores[i]->UC[j][k] = isochores[0]->UC[j][k];
      
	  for(k=0; k< isochores[i]->ne[j]; k++)
	    isochores[i]->DE[j][k] = isochores[0]->DE[j][k];

	  isochores[i]->block[j]  = isochores[0]->block[j];
	}
    }
}

char* RequestMemorySequence(long L)
{
  char* s;

  if ((s = (char*) calloc(L,sizeof(char))) == NULL)
    printError("Not enough space to hold ADN sequence"); 

  return(s);
}

packSites* RequestMemorySites()
{
  packSites* allSites;

  /* Allocating memory for sites */
  if ((allSites = 
       (struct s_packSites *) malloc(sizeof(struct s_packSites))) == NULL)
    printError("Not enough space to hold pack of sites");  
  
  /* Starts */
  if ((allSites->StartCodons = 
       (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough space to hold start codons");
  
  /* AcceptorSites */
  if ((allSites->AcceptorSites = 
       (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough space to hold acceptor sites");
	      	      
  /* DonorSites */
  if ((allSites->DonorSites = 
       (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough space to hold donor sites");
	      
  /* Stop Codons */
  if ((allSites->StopCodons = 
       (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough space to hold stop codons");
  /* end memory-sites allocation */

  return(allSites);
}

packExons* RequestMemoryExons()
{
  packExons* allExons;

  /* Allocating memory for exons */
  if ((allExons = 
       (struct s_packExons *) malloc(sizeof(struct s_packExons))) == NULL)
    printError("Not enough space to hold pack of exons");  

  /* InitialExons */
  if ((allExons->InitialExons = 
       (exonGFF*) calloc(NUMEXONS, sizeof(exonGFF))) == NULL)
    printError("Not enough space to hold initial exons");
      
  /* InternalExons */
  if ((allExons->InternalExons = 
       (exonGFF*) calloc(NUMEXONS, sizeof(exonGFF))) == NULL)
    printError("Not enough space to hold internal exons");
          
  /* TerminalExons */
  if ((allExons->TerminalExons = 
       (exonGFF*) calloc(NUMEXONS, sizeof(exonGFF))) == NULL)
    printError("Not enough space to hold terminal exons");

  /* TerminalExons */
  if ((allExons->Singles = 
       (exonGFF*) calloc(NUMEXONS, sizeof(exonGFF))) == NULL)
    printError("Not enough space to hold single genes");
  /* end memory-exons allocation */

  return(allExons);
}

exonGFF* RequestMemorySortExons()
{
  exonGFF *exons;
  /* Sorting Exons */
  if ((exons =
       (exonGFF*) calloc(NUMEXONS * 4, sizeof(exonGFF)))  == NULL)
    printError("Not enough space to hold all exons(sorting)");

  return(exons);
}

packEvidence* RequestMemoryEvidence()
{
  packEvidence* p;

  /* Allocating memory for structure */
  if ((p = 
       (struct s_packEvidence *) malloc(sizeof(struct s_packEvidence))) == NULL)
    printError("Not enough space to hold pack of evidence exons");    

  /* Evidence sites */
  if ((p->vSites = 
       (struct s_site *) calloc(NUMSEVIDENCES, sizeof(struct s_site))) == NULL)
    printError("Not enough space to hold evidence sites");

  /* Evidence exons */
  if ((p->vExons =
       (exonGFF*) calloc(NUMEEVIDENCES, sizeof(exonGFF)))  == NULL)
    printError("Not enough space to hold evidence exons");


  /* Reset counters */
  p->nvSites = 0;
  p->nvExons = 0;
  p->i1vExons = 0;
  p->i2vExons = 0;

  return(p);
}

packSR* RequestMemorySimilarityRegions()
{
  packSR* p;
  int i;

  /* Allocating memory for structure */
  if ((p = 
       (struct s_packSR *) malloc(sizeof(struct s_packSR))) == NULL)
    printError("Not enough space to hold pack of similarity regions");    

  /* Similarity regions */
  if ((p->sRegions = 
       (SR **) calloc(STRANDS*FRAMES, sizeof(SR*))) == NULL)
    printError("Not enough space to hold similarity regions array");

  for(i=0;i<STRANDS*FRAMES;i++)
    if ((p->sRegions[i] = 
	 (SR *) calloc(MAXSR, sizeof(SR))) == NULL)
      printError("Not enough space to hold similarity regions array");  
  
  /* Counters */
  if ((p->nRegions =
       (int*) calloc(STRANDS*FRAMES, sizeof(int)))  == NULL)
    printError("Not enough space to hold similarity regions counters");
  if ((p->iRegions =
       (int*) calloc(STRANDS*FRAMES, sizeof(int)))  == NULL)
    printError("Not enough space to hold similarity regions counters");


  /* Reset counters */
  for(i=0;i<STRANDS*FRAMES;i++)
    {
      p->nRegions[i] = 0;
      p->iRegions[i] = 0;
    }

  return(p);
}

gparam* RequestMemoryParams()
{
  gparam* gp;
  long OligoDim;

  gp =(gparam *) malloc(sizeof(gparam));

  /* Profiles */
  if ((gp->StartProfile = (profile *) malloc(sizeof(profile))) == NULL)  
    printError("Not enough space to hold profiles");
  if ((gp->AcceptorProfile = (profile *) malloc(sizeof(profile))) == NULL)  
    printError("Not enough space to hold profiles");
  if ((gp->DonorProfile = (profile *) malloc(sizeof(profile))) == NULL)  
    printError("Not enough space to hold profiles");
  if ((gp->StopProfile = (profile *) malloc(sizeof(profile))) == NULL)  
    printError("Not enough space to hold profiles");

  /* Oligo arrays (Markov models) */
  OligoDim = (int)pow((float)4,(float)OLIGOLENGTH);
  
  if ((gp->OligoLogsIni[0]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough space to hold Markov pentamers");

  if ((gp->OligoLogsIni[1]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough space to hold Markov pentamers");

  if ((gp->OligoLogsIni[2]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough space to hold Markov pentamers");

  if ((gp->OligoLogsTran[0]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough space to hold Markov hexamers");

  if ((gp->OligoLogsTran[1]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough space to hold Markov hexamers");

  if ((gp->OligoLogsTran[2]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough space to hold Markov hexamers");

  /* Score parameters */
  if ((gp->Initial = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough space to hold scoring parameters");

  if ((gp->Internal = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough space to hold scoring parameters");
 
  if ((gp->Terminal = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough space to hold scoring parameters");

  if ((gp->Single = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough space to hold scoring parameters");

  return(gp);
}

gparam ** RequestMemoryIsochoresParams()
{
  gparam** isochores;
  int i;

  /* Allocating the array of isochores */
  if ((isochores = (gparam **) calloc(MAXISOCHORES, sizeof(gparam *))) == NULL)
    printError("Not enough space to hold isochores array");

  /* Allocating each position */
  for(i=0; i<MAXISOCHORES; i++)
    isochores[i] = (gparam *) RequestMemoryParams();
  
  /* Allocating space for global parameters */
  if ((isochores[0]->D = (dict *)malloc(sizeof(dict))) == NULL)
    printError("Not enough space to hold ExonType-Dictionary");

  if ((isochores[0]->nc = (int *) calloc(MAXENTRY, sizeof(int))) == NULL)
    printError("Not enough space to hold nc-array");

  if ((isochores[0]->ne = (int *) calloc(MAXENTRY, sizeof(int))) == NULL)
    printError("Not enough space to hold ne-array");

  if ((isochores[0]->md = (long *) calloc(MAXENTRY, sizeof(long))) == NULL)
    printError("Not enough space to hold minDist-array");

  if ((isochores[0]->Md = (long *) calloc(MAXENTRY, sizeof(long))) == NULL)
    printError("Not enough space to hold max-Dist-array");

  return(isochores);
}


void RequestMemoryProfile(profile* p)
{
  int i;
  
  /* Transition probabilities */
  p->dimensionTrans = (int)pow((float)4,(float)(p->order+1));   
  for(i=0; i<PROFILEDIM; i++)
    if ((p->transitionValues[i] = 
	 (float *) calloc(p->dimensionTrans, sizeof(float))) == NULL)
      printError("Not enough space to hold this profile");
}

packGenes* RequestMemoryGenes()
{
  packGenes* pg;
  int aux, aux2;

  /* Allocating memory for packGenes */
  if ((pg = 
       (struct s_packGenes *) malloc(sizeof(struct s_packGenes))) == NULL)
    printError("Not enough space to hold pack of genes");  

  /* Allocating memory space for Ghost Exon */
  if (( pg->Ghost = (exonGFF *) malloc(sizeof(exonGFF))) == NULL)
    printError("Not enough space to hold Ghost Exon");

  pg->Ghost->Acceptor = (site*)calloc(1,sizeof(site));
  pg->Ghost->Donor = (site*)calloc(1,sizeof(site));

  /* Mark exon as Ghost Exon */
  pg->Ghost -> Strand = '*';
  pg->Ghost -> GeneScore = 0.0;
  pg->Ghost -> Donor -> Position = 0;
  pg->Ghost -> offset1 = 0;
  pg->Ghost -> offset2 = 0;
  pg->GOptim = pg->Ghost;
  
  if ((pg->Ga = (exonGFF* **)calloc(MAXENTRY, sizeof(exonGFF* *))) == NULL)
    printError("Not enough space to hold Ga-Exons");

  /* Initialize Ga-exons: All look at the Ghost exon */
  for(aux=0; aux<MAXENTRY; aux++) 
    {
      if ((pg->Ga[aux] = (exonGFF* *)calloc(FRAMES, sizeof(exonGFF*))) == NULL)
	printError("Not enough space to hold frame of Ga-Exons");

      for(aux2=0; aux2 < FRAMES; aux2++)
      	pg->Ga[aux][aux2] = pg->Ghost;
    }

  /* Alloc memory space for the set of arrays */
  /* Memory for the array of sorting functions */
  if ((pg->d = (exonGFF* **)calloc(MAXENTRY, sizeof(exonGFF* *))) == NULL)
    printError("Not enough space to hold set of arrays");
  
  /* Memory for each sorting function */
  for(aux=0; aux < MAXENTRY; aux++) 
    if ((pg->d[aux] = (exonGFF* *)calloc(2 * NUMEXONS, sizeof(exonGFF*))) == NULL)
      printError("Not enough space to hold sorted array of pointer to exons");
  
  if ((pg->km = (long *)calloc(MAXENTRY, sizeof(long))) == NULL)
    printError("Not enough space to hold counters of sorting functions");

  if ((pg->je = (long *)calloc(MAXENTRY, sizeof(long))) == NULL)
    printError("Not enough space to hold j-counter of each class");

  return(pg);
}

packDump* RequestMemoryDumpster()
{
  packDump* d;

  /* Allocating memory for dumpster */
  if ((d = 
       (struct s_packDump *) malloc(sizeof(struct s_packDump))) == NULL)
    printError("Not enough space to hold dumpster");  

  /* Dumpster Sites(bestgenes) */
  if ((d->dumpSites = 
       (struct s_site *) calloc(MAXBACKUPSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough space to hold dumpster sites");

  /* Dumpster exons(bestgenes) */
  if ((d->dumpExons = 
       (exonGFF*) calloc(MAXBACKUPEXONS, sizeof(struct s_exonGFF))) == NULL)
    printError("Not enough space to hold dumpster exons");  

  /* Dumpster hash */
  if ((d->h = 
       (dumpHash*) malloc(sizeof(struct s_dumpHash))) == NULL)
    printError("Not enough space to hold dumpster hash table");  
  
  if ((d->h->T = 
       (dumpNode**) calloc(MAXBACKUPEXONS,sizeof(dumpNode*))) == NULL)
    printError("Not enough space to hold dumpster hash table");   

  d->ndumpSites = 0;
  d->ndumpExons = 0;
  resetDumpHash(d->h);
  
  return(d);
}

dict* RequestMemoryAaDictionary()
{
   dict* dAA;

  /* 1. Allocating memory for aminoacid dictionary */
  if ((dAA = (dict *)malloc(sizeof(dict))) == NULL)
    printError("Not enough space to hold Aa-Dictionary");

  /* 2. Reset dictionary */
  resetDict(dAA);

  /* 3. Filling it with the aminoacids and the codons */
  setAADict(dAA,"GCA",'A');
  setAADict(dAA,"GCC",'A');
  setAADict(dAA,"GCG",'A');
  setAADict(dAA,"GCT",'A');
  setAADict(dAA,"TGC",'C');
  setAADict(dAA,"TGT",'C');
  setAADict(dAA,"GAC",'D');
  setAADict(dAA,"GAT",'D');
  setAADict(dAA,"GAA",'E');
  setAADict(dAA,"GAG",'E');
  setAADict(dAA,"TTC",'F');
  setAADict(dAA,"TTT",'F');
  setAADict(dAA,"GGA",'G');
  setAADict(dAA,"GGC",'G');
  setAADict(dAA,"GGG",'G');
  setAADict(dAA,"GGT",'G');
  setAADict(dAA,"CAC",'H');
  setAADict(dAA,"CAT",'H');
  setAADict(dAA,"ATA",'I');
  setAADict(dAA,"ATC",'I');
  setAADict(dAA,"ATT",'I');
  setAADict(dAA,"AAA",'K');
  setAADict(dAA,"AAG",'K');
  setAADict(dAA,"TTA",'L');
  setAADict(dAA,"TTG",'L');
  setAADict(dAA,"CTA",'L');
  setAADict(dAA,"CTC",'L');
  setAADict(dAA,"CTG",'L');
  setAADict(dAA,"CTT",'L');
  setAADict(dAA,"ATG",'M');
  setAADict(dAA,"AAC",'N');
  setAADict(dAA,"AAT",'N');
  setAADict(dAA,"CCA",'P');
  setAADict(dAA,"CCC",'P');
  setAADict(dAA,"CCG",'P');
  setAADict(dAA,"CCT",'P');
  setAADict(dAA,"CAA",'Q');
  setAADict(dAA,"CAG",'Q');
  setAADict(dAA,"AGA",'R');
  setAADict(dAA,"AGG",'R');
  setAADict(dAA,"CGA",'R');
  setAADict(dAA,"CGC",'R');
  setAADict(dAA,"CGG",'R');
  setAADict(dAA,"CGT",'R');
  setAADict(dAA,"AGC",'S');
  setAADict(dAA,"AGT",'S');
  setAADict(dAA,"TCA",'S');
  setAADict(dAA,"TCC",'S');
  setAADict(dAA,"TCG",'S');
  setAADict(dAA,"TCT",'S');
  setAADict(dAA,"ACA",'T');
  setAADict(dAA,"ACC",'T');
  setAADict(dAA,"ACG",'T');
  setAADict(dAA,"ACT",'T');
  setAADict(dAA,"GTA",'V');
  setAADict(dAA,"GTC",'V');
  setAADict(dAA,"GTG",'V');
  setAADict(dAA,"GTT",'V');
  setAADict(dAA,"TGG",'W');
  setAADict(dAA,"TAC",'Y');
  setAADict(dAA,"TAT",'Y');
  setAADict(dAA,"TAA",'*');
  setAADict(dAA,"TAG",'*');
  setAADict(dAA,"TGA",'*');

  return(dAA);
}


