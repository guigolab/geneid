/*************************************************************************
*                                                                        *
*   Module: geneid                                                       *
*                                                                        *
*   Main program. Management of the actions of geneid.                   *
*                                                                        *
*   This file is part of the geneid Distribution                         *
*                                                                        *
*     Copyright (C) 2000 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           * 
*     with contributions from:                                           *
*                          Moises  BURSET ALVAREDA                       *
*                          Xavier  MESSEGUER PEYPOCH                     *
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

/* $Id: geneid.c,v 1.1 2000-07-18 18:03:35 eblanco Exp $ */

#include "geneid.h"

/* Setup Flags of geneid */
int   SFP=0, SDP=0, SAP=0, STP=0,
      EFP=0, EIP=0, ETP=0, EXP=0, ESP=0,
      VRB=0, FWD=1, RVS=1,
      GENAMIC = 1, GENEID = 1, 
      GFF = 0, X10 = 0, 
      EVD = 0, SRP = 0;

float EW = NOVALUE;

/* Accounting time and results */
account *m;

long analizeFile(char* SequenceFile)
{
  struct stat *buf;
  long size;

  buf = (struct stat*)malloc(sizeof(struct stat));

  if (stat(SequenceFile,buf) != 0)
    printError("Impossible to read sequence file");

  size = (long)buf->st_size;
  
  if (size == 0)
    printError("Empty sequence file!");

  return(size);
}

int SelectIsochore(float percent, gparam** isochores)
{
  int i;
  int stop;

  i = 0;

  percent = PERCENT * percent;

  stop = (percent >= isochores[i]->leftValue && 
	  percent <= isochores[i]->rightValue);

  while(!stop)
    {
      i++;
      stop = (percent >= isochores[i]->leftValue && 
	      percent <= isochores[i]->rightValue);
    }

  return(i);
}

/* Management of splice sites prediction and exon building */
void  PredictExons(char *Sequence, long LengthSequence,
		   packSites* allSites, packExons* allExons,
		   long l1, long l2,
		   int Strand, packSR* sr, 
		   gparam* gp,
		   gparam** isochores, int nIsochores)
{
  char mess[MAXSTRING];
  
  long l1a, l1b,
       l2a, l2b,
       l1c, l2c;

  long cutPoint;
  
  /* 0. Define boundaries of splice site prediction
        according to current split positions and strand selected */
  if (Strand == FORWARD)
    {
      /* Forward direction */
      /* Starts and Acc boundaries */
      l1a = l1;
      l2a = (l2 == LengthSequence - 1)? l2 : l2 - OVERLAP;

      /* Donors Boundary */
      l1b = l1;
      l2b = l2;

      /* Stops Boundary */
      l1c = l1;
      l2c = l2;

      /* Exons can be built with Stops before cutPoint */
      cutPoint = l1;
    }
  else
    {
      /* Reverse direction */
      /* Starts and Acc Boundaries */
      l1a = l1;
      l2a = l2;

      /* Donor Boundary */
      l1b = (l1 == 0)? l1: l1 + OVERLAP;
      l2b = l2;

      /* Stops boundary */
      l1c = l1;
      l2c = l2;

      /* Exons can be built with Stops before cutPoint */
      cutPoint = (l1 == 0)? l1 : l1 + OVERLAP;
    }

  /* 1. Predicting splice sites of current split of DNA sequence */ 
  printMess ("Computing sites ...");

  allSites->nStartCodons =
    GetSitesWithProfile(Sequence,gp->StartProfile,allSites->StartCodons,l1a,l2a);
  sprintf(mess, "Start Codons \t\t%8ld", allSites->nStartCodons);
  printRes(mess);

  allSites->nAcceptorSites =
    GetSitesWithProfile(Sequence,gp->AcceptorProfile,allSites->AcceptorSites,l1a,l2a);
  sprintf(mess, "Acceptor Sites \t\t%8ld", allSites->nAcceptorSites);
  printRes(mess);

  allSites->nDonorSites =
    GetSitesWithProfile(Sequence,gp->DonorProfile,allSites->DonorSites,l1b,l2b);
  sprintf (mess,"Donor Sites \t\t%8ld", allSites->nDonorSites);
  printRes(mess);

  allSites->nStopCodons =
  GetStopCodons(Sequence,gp->StopProfile, allSites->StopCodons,l1c,l2c);
  sprintf (mess,"Stop Codons \t\t%8ld", allSites->nStopCodons);
  printRes(mess);
  
  /* Total number of predicted splice sites in this strand */
  allSites->nSites = allSites->nStartCodons +
      allSites->nAcceptorSites +
      allSites->nDonorSites +
      allSites->nStopCodons;

  sprintf(mess,"---------\t\t%8ld", allSites->nSites);
  printRes(mess);

  /* 2. Building exons with splice sites predicted before */ 
  printMess ("Computing exons ...");   
  
  allExons->nInitialExons =
    BuildInitialExons(allSites->StartCodons,allSites->nStartCodons,
		      allSites->DonorSites,allSites->nDonorSites,
		      allSites->StopCodons,allSites->nStopCodons,
		      gp->MaxDonors,allExons->InitialExons);
  sprintf(mess,"Initial Exons \t\t%8ld", allExons->nInitialExons);
  printRes(mess); 
  
  allExons->nInternalExons =
    BuildInternalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
		       allSites->DonorSites,allSites->nDonorSites,
		       allSites->StopCodons,allSites->nStopCodons,
		       gp->MaxDonors,allExons->InternalExons);
  sprintf(mess,"Internal Exons \t\t%8ld", allExons->nInternalExons);
  printRes(mess); 
  
  allExons->nTerminalExons =
    BuildTerminalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
             allSites->StopCodons,allSites->nStopCodons,
             LengthSequence,allExons->TerminalExons,
             cutPoint);
  sprintf(mess,"Terminal Exons \t\t%8ld", allExons->nTerminalExons);
  printRes(mess); 

  allExons->nSingles =
    BuildSingles(allSites->StartCodons,allSites->nStartCodons,
		 allSites->StopCodons,allSites->nStopCodons,
		 cutPoint, allExons->Singles);
  sprintf(mess,"Single genes \t\t%8ld", allExons->nSingles);
  printRes(mess); 

  /* 3. Scoring and Filtering Exons */
  
  /* Scoring and filtering built exons */
  ScoreExons(Sequence, allExons, 
	     l1, l2, Strand, sr,
	     isochores,nIsochores); 
  
  /* Total number of built exons in this strand */
  allExons->nExons = allExons->nInitialExons + 
         allExons->nInternalExons +
         allExons->nTerminalExons +
         allExons->nSingles;

  sprintf(mess,"---------\t\t%8ld", allExons->nExons);
  printRes(mess); 
}

/*****************************************************************************/

int main (int argc, char *argv[])
{
  /* DNA-sequence */
  FILE* seqfile;
  char* Sequence;
  char* RSequence;
  long  LengthSequence;

  /* Boundaries of current split */
  long l1,l2;

  /* Forward Data Structures */
  packSites* allSites;
  packExons* allExons;

  /* Reverse Data Structures */
  packSites* allSites_r;
  packExons* allExons_r;

  /* Sort by acceptor Exons */
  exonGFF* exons; 
  long nExons;

  /* Evidence exons */
  packEvidence* evidence;
  long currVExons=0;

  /* Similarity to protein regions */
  packSR* sr;

  /* Geneid prediction parameters */
  gparam* gp;
  gparam** isochores;
 
  /* Input Files */
  char   SequenceFile[FILENAMELENGTH],
         ExonsFile[FILENAMELENGTH],
         SRFile[FILENAMELENGTH],
         ParamFile[FILENAMELENGTH]="";
  
  /* Name of current sequence */ 
  char Locus[LOCUSLENGTH];
  char nextLocus[LOCUSLENGTH];
  
  char mess[MAXSTRING];

  int lastSplit; 

  /* Best partial predicted genes */
  packGenes* genes;   

  /* Dumpster for back-up operations between splits */
  packDump* dumpster;

  /* Amino-acid dictionary */
  dict* dAA;

  /* Measure of C+G content */
  float percent;
  int currentIsochore;
  int nIsochores;

  int reading;

  /* 0. Read setup options */
  readargv(argc,argv,ParamFile,SequenceFile,ExonsFile,SRFile);

  LengthSequence = analizeFile(SequenceFile);
  sprintf(mess,"ADN sequence file size = %ld bytes",LengthSequence);
  printMess(mess);

  /* 1. Alloc main memory structures */
  printMess("Request Memory to System");
  Sequence   = (char*)         RequestMemorySequence(LengthSequence);
  RSequence  = (char*)         RequestMemorySequence(LengthSequence);
  allSites   = (packSites*)    RequestMemorySites();
  allSites_r = (packSites*)    RequestMemorySites();
  allExons   = (packExons*)    RequestMemoryExons();
  allExons_r = (packExons*)    RequestMemoryExons();
  exons      = (exonGFF*)      RequestMemorySortExons();
  evidence   = (packEvidence*) RequestMemoryEvidence();
  genes      = (packGenes*)    RequestMemoryGenes();
  dumpster   = (packDump*)     RequestMemoryDumpster();
  dAA        = (dict*)         RequestMemoryAaDictionary();
  sr         = (packSR*)       RequestMemorySimilarityRegions();
  isochores  = (gparam**)      RequestMemoryIsochoresParams();

  /* 2. Read prediction parameters file */
  printMess("Reading parameters..."); 
  nIsochores = readparam(ParamFile, isochores);

  /* initializing stats... */
  m = (account*)InitAcc();  

  /* 3. Predict Gene or Predict Exons? */
  if (GENEID)
    {
      /* a. Predicting splice sites, exons and genes of DNA sequences */
      
      /* a.0. Reading ADN sequence */
      if ((seqfile = fopen(SequenceFile, "rb"))==NULL) 
        printError("The Sequence file can not be open for read");

      /* reading first sequence (in Fasta format) */
      reading = IniReadSequence(seqfile,Locus);
      while (reading!=EOF)
	{
	  printMess("Reading ADN sequence");
	  reading = ReadSequence(seqfile, Sequence, nextLocus);
	  
	  /* a.1. Processing sequence: Reverse and complement */
	  printMess("Processing ADN sequence");
	  LengthSequence = FetchSequence(Sequence, RSequence);

	  if (!LengthSequence)
	    printError("Sequence empty");
 	 
	  /* Header Output */
	  OutputHeader(Locus, LengthSequence);
 	  
	  /* a.2. Reading evidence exons */
	  if (EVD)
	    { 
	      printMess("Reading evidence exons");
	      evidence->nvExons = ReadExonsGFF(ExonsFile, evidence, isochores[0]->D); 
	      sprintf(mess,"%ld evidence exons read\n", evidence->nvExons); 
	      printMess(mess);
	    }

	  /* a.2. Reading Similarity to protein regions */
	  if (SRP)
	    {
	      printMess("Reading Similarity to protein regions");
	      sprintf(mess,"%ld Similarity to protein regions read\n", 
		      ReadSR(SRFile, sr, LengthSequence)); 
	      printMess(mess); 
	    }
  
	  /* a.3. Split sequence in several parts of the same length */
	  /* l1 is the left boundary and l2 is the right */
	  l1 = 0;
	  l2 = MIN(LENGTHSi-1,LengthSequence-1);
	  lastSplit = (l2 == LengthSequence-1);
	  while((l1 < (LengthSequence - OVERLAP))|| (l1 == 0))
	    {
	      /* Measure of C+G content of this region */
	      percent = MeasureSequence(l1,l2,Sequence);
	      sprintf(mess,"C+G-Percent is %f\n", percent);
	      printMess(mess); 
	      
	      currentIsochore = SelectIsochore(percent,isochores);
	      gp = isochores[currentIsochore];
	      sprintf(mess,"Selecting isochore %d\n", currentIsochore+COFFSET);
	      printMess(mess);

	      /* a.4. Prediction of sites and exons building */
	      if (FWD) 
		{
		  /* forward sense predictions */ 
		  sprintf(mess,"Running FWD  %s: %ld - %ld", Locus,l1,l2);
		  printMess(mess);
		  PredictExons(Sequence, LengthSequence, 
			       allSites, allExons, l1, l2, FORWARD, sr, gp,
			       isochores,nIsochores);
		}      
	      if (RVS) 
		{
		  /* reverse sense predictions */ 
		  sprintf(mess,"Running Reverse  %s: %ld - %ld(%ld - %ld)", 
			  Locus, LengthSequence-1 -l2, 
			  LengthSequence-1 -l1,l1,l2);         
		  printMess(mess);
		  PredictExons(RSequence, LengthSequence,
			       allSites_r, allExons_r, 
			       LengthSequence-1 - l2, 
			       LengthSequence-1 - l1, REVERSE, sr, gp,
			       isochores,nIsochores);
		  
		  /* computing sites positions according to forward sense positions */
		  RecomputePositions(allSites_r, LengthSequence);
		  
		  /* exchange acc and donor sites because of Acc < Don needed */
		  SwitchPositions(allExons_r);
		}
	      
	      /* a.5. Sort (forward & reverse) exons by acc position */
	      if (EVD)
		{
		  /* Searching evidence exons in this split */
		  printMess("Searching evidence exons");
		  currVExons = SearchEvidenceExons(evidence, 
				      (lastSplit)?l2:l2-OVERLAP);
		}
	      
	      nExons = allExons->nExons + allExons_r->nExons + currVExons;
	      sprintf(mess,"Sorting %ld Exons by Acceptor\n", nExons);  
	      printMess(mess);
  	      
	      /* Merge predicted exons with some evidence exons */
	      SortExons(allExons, allExons_r, evidence, exons, l1, l2);
	      
	      if (EVD)
		SwitchCounters(evidence);

	      /* a.6. Printing current split predictions (sites and exons) */
	      Output(allSites, allSites_r, allExons, allExons_r, 
		     exons, nExons, Locus, l1, l2, Sequence, gp, dAA); 
	      
	      /* recompute stats about splice sites and exons */
	      updateTotals(m,allSites,allSites_r,allExons,allExons_r);
	      
	      /* a.7. Calling genamic for assembling the best gene */ 
	      if (GENAMIC && nExons)
		{
		  genamic(exons, nExons, genes, gp);
		  printMess("Ending GenAmic...");
		  
		  /* clean hash table of exons */
		  cleanDumpHash(dumpster->h);

		  /* a.8. Backup operations of genes to next split */
		  if (!lastSplit)
		    {
		      /* backup of unused genes */
		      printMess("Backup of d-genes");
		      BackupArrayD(genes, l2 - OVERLAP, gp, dumpster);
		      
		      /* back-up best partial genes */
		      printMess("Backup of partial genes");
		      BackupGenes(genes, gp->nclass, dumpster);
		    }
		}
	      
	      /* a.9. Compute boundaries of next split of sequence */
	      l1 += LENGTHSi - OVERLAP;
	      l2 = MIN(l1 + LENGTHSi -1, LengthSequence-1);
	      lastSplit = (l2 == LengthSequence-1);
	    }/* End of current sequence */
	  
	  if (GENAMIC)
	    {	      
	      /* a.10. Printing gene predictions */
	      OutputGene(genes,m->totalExons, Locus, Sequence, gp, dAA);
	      
	      /* reset backup genes data structures to the next DNA sequence */
	      printMess("Cleanning genes");
	      cleanGenes(genes,gp->nclass,dumpster);
	    }
	  
	  /* showing global stats about last sequence predicted */
	  OutputStats(Locus, evidence->nvExons);
	  cleanAcc(m);
	  strcpy(Locus,nextLocus);
	} /* endwhile(reading) */
      
    } /*endifgeneid*/
  else
    {
      /* b. Predicting gens of the exonsGFF using blocks */
      /* It's necessary to increase NUMSITES and NUMEXONS */

      /* b.0. Reading ADN sequence */
      /* open the Sequence File */
      if ((seqfile = fopen(SequenceFile, "rb"))==NULL) 
        printError("The Sequence file can not be open for read");
      
      /* read the sequence (because of translation to proteins) */
      printMess("Reading ADN sequence");
      reading = IniReadSequence(seqfile,Locus);
      if (reading != EOF)
	{
	  reading = ReadSequence(seqfile, Sequence, nextLocus);
	  LengthSequence = FetchSequence(Sequence, RSequence);
	}

      /* Header Output */
      OutputHeader(Locus, LengthSequence);
      
      /* b.1. Reading exons in GFF format */
      printMess("Reading exonsGFF from file");
      nExons = ReadExonsGFF(ExonsFile, evidence, isochores[0]->D); 
      
      sprintf(mess,"%ld exons GFF stored\n", nExons); 
      printMess(mess); 
      
      printMess("Executing GenAmic..."); 
      if (nExons) 
	{ 
	  /* b.2. Calling genamic for assembling the best gene */
	  genamic(evidence->vExons, nExons, genes, isochores[0]);  
	  printMess("Ending GenAmic..."); 
	}
      
      /* b.3. Printing gene predictions */
      OutputGene(genes, nExons, Locus, Sequence, isochores[0], dAA);
    } /* endGenamic */

  /* 4. The End */
  OutputTime(); 
  exit(0);
  return(0);
}
