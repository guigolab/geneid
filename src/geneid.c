/*************************************************************************
*                                                                        *
*   Module: geneid                                                       *
*                                                                        *
*   geneid main program                                                  *
*                                                                        *
*   This file is part of the geneid 1.1 distribution                     *
*                                                                        *
*     Copyright (C) 2001 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           * 
*     with contributions from:                                           *
*                          Moises  BURSET ALVAREDA                       *
*                          Genis   PARRA FARRE                           *
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

/* $Id: geneid.c,v 1.10 2003-02-26 10:48:14 eblanco Exp $ */

#include "geneid.h"

/* geneid setup flags */
int
  /* sites to print */
  SFP=0, SDP=0, SAP=0, STP=0,
  /* exons to print */
  EFP=0, EIP=0, ETP=0, EXP=0, ESP=0, EOP = 0,
  /* Partial or full prediction engine */
  GENAMIC = 1, GENEID = 1, 
  /* Only forward or reverse prediction engine */
  FWD=1, RVS=1,
  /* switch ORF prediction on */
  scanORF = 0,
  /* Input annotations or homology to protein information */
  EVD = 0, SRP = 0,
  /* Output formats */
  GFF = 0, X10 = 0, XML = 0, cDNA = 0,
  /* Verbose flag (memory/processing information) */
  BEG=0, VRB=0;

/* Increase/decrease exon weight value (exon score) */
float EW = NOVALUE;

/* Generic maximum values: sites, exons and backup elements */
long NUMSITES,NUMEXONS,MAXBACKUPSITES,MAXBACKUPEXONS;

/* Accounting time and results */
account *m;

/************************************************************************
                            geneid MAIN program
************************************************************************/

int main (int argc, char *argv[])
{
  /* DNA sequence data structures */
  FILE* seqfile;
  char* Sequence;
  char* RSequence;
  long  LengthSequence;
  
  /* Current split ends */
  long l1,l2;
  
  /* Forward semse data structures */
  packSites* allSites;
  packExons* allExons;
  
  /* Reverse sense data structures */
  packSites* allSites_r;
  packExons* allExons_r;
  
  /* Table to sort predicted exons by acceptor */
  exonGFF* exons; 
  long nExons;
  
  /* Evidences (annotations) input by user */
  packEvidence* evidence;
  
  /* Similarity to protein regions input by user */
  packSR* sr;
  
  /* Best partial predicted genes */
  packGenes* genes;   
  
  /* Dumpster for backup operations between splits */
  packDump* dumpster;
  
  /* Amino acid dictionary (genetic code) */
  dict* dAA;
  
  /* geneid prediction parameters: data structures */
  gparam* gp;
  gparam** isochores;
  
  /* Input Filenames */
  char   SequenceFile[FILENAMELENGTH],
	ExonsFile[FILENAMELENGTH],
	SRFile[FILENAMELENGTH],
	ParamFile[FILENAMELENGTH]="";
  
  /* Locus sequence name */ 
  char Locus[LOCUSLENGTH];
  char nextLocus[LOCUSLENGTH];
  
  /* Measure of C+G content to select the isochore */
  packGC* GCInfo;
  packGC* GCInfo_r;
  int inigc, endgc;
  float percentGC;
  int currentIsochore;
  int nIsochores;
  
  int lastSplit;
  char mess[MAXSTRING];
  int reading;
  
  
  /** 0. Starting and reading options, parameters and sequence... **/
  
  /* 0.a. Previous checkpoint about length in splits and overlapping */
  if (LENGTHSi <= OVERLAP)
    printError("LENGTHSi must be greater than OVERLAP parameter (geneid.h)");
  
  /* 0.b. Initializing stats and time counters */
  m = (account*)InitAcc();  
  
  /* 0.c. Read setup options */
  readargv(argc,argv,ParamFile,SequenceFile,ExonsFile,SRFile);
  printRes("\n\n\t\t\t** Executing geneid 1.1 2001 geneid@imim.es **\n\n");
  
  /* 0.d. Prediction of DNA sequence length to request memory */
  LengthSequence = analizeFile(SequenceFile);
  sprintf(mess,"DNA sequence file size = %ld bytes",LengthSequence);
  printMess(mess);
  
  /* 0.e. Computing ratios for every type of signal and exons */
  printMess("Computing Ratios");
  SetRatios(&NUMSITES,
            &NUMEXONS,
            &MAXBACKUPSITES,
            &MAXBACKUPEXONS,
            LengthSequence);
  
  /** 1. Allocating main geneid data structures **/
  printMess("Request Memory to Operating System");
  
  /* 1.a. Mandatory geneid data structures */
  Sequence   = (char*)         RequestMemorySequence(LengthSequence);
  RSequence  = (char*)         RequestMemorySequence(LengthSequence);
  allSites   = (packSites*)    RequestMemorySites();
  allSites_r = (packSites*)    RequestMemorySites();
  allExons   = (packExons*)    RequestMemoryExons();
  allExons_r = (packExons*)    RequestMemoryExons();
  exons      = (exonGFF*)      RequestMemorySortExons();
  genes      = (packGenes*)    RequestMemoryGenes();
  isochores  = (gparam**)      RequestMemoryIsochoresParams(); 
  GCInfo     = (packGC*)       RequestMemoryGC();
  GCInfo_r   = (packGC*)       RequestMemoryGC();
  dAA        = (dict*)         RequestMemoryAaDictionary();
  
  /* 1.b. Optional data structures according to the selected options */
  if (EVD || !GENEID)
    evidence   = (packEvidence*) RequestMemoryEvidence();
  else
	evidence = NULL;

  if (SRP)
    sr         = (packSR*)       RequestMemorySimilarityRegions();
  else
	sr = NULL;

  /* 1.c. Backup information might be necessary between splits */
  if (LengthSequence > LENGTHSi)
    dumpster   = (packDump*)     RequestMemoryDumpster();
  else
	dumpster = NULL;

  /** 2. Reading statistical model parameters file **/
  printMess("Reading parameters..."); 
  nIsochores = readparam(ParamFile, isochores);
  
  /** 3. Starting processing: complete or partial prediction ? **/
  if (GENEID)
    {
      /* A. Predicting signals, exons and genes in DNA sequences */
      /* A.1. Reading evidences (annotations) */
	  if (EVD)
		{ 
		  printMess("Reading evidences (annotations)");
		  evidence->nvExons =
			ReadExonsGFF(ExonsFile, evidence, isochores[0]->D);
		  sprintf(mess,"%ld evidences read\n", evidence->nvExons); 
		  printMess(mess);
		}
	
	  /* MOVING THIS TO THE NEXT LOOP */
  
      /** A.3. Input DNA sequences (perhaps more than one) **/
      if ((seqfile = fopen(SequenceFile, "rb"))==NULL) 
        printError("The input sequence file can not be opened to read");
	  
      /* reading the locusname of sequence (in Fasta format) */
      reading = IniReadSequence(seqfile,Locus);
	  
      while (reading != EOF)
		{		  
          printMess("Loading DNA sequence");
		  reading = ReadSequence(seqfile, Sequence, nextLocus);
		  		  
          /* A.4. Prepare sequence to work on */
          printMess("Processing DNA sequence");
          LengthSequence = FetchSequence(Sequence, RSequence);
		
		  /* A.2. Reading homology to protein information */
		  if (SRP)
			{
			  printMess("Reading similarity regions (homology information)");
			  sr->nTotalRegions = ReadSR(SRFile, sr, LengthSequence);
			  sprintf(mess,"%d similarity regions (SR) read\n",sr->nTotalRegions);
			  printMess(mess); 
			}
  
          /* Output header information */
          OutputHeader(Locus, LengthSequence);

		  /* A.5. Processing sequence into several fragments if required */
		  /* l1 is the left end and l2 is the right end in Sequence */
		  l1 = 0;
		  l2 = MIN(LENGTHSi-1,LengthSequence-1);
		  lastSplit = (l2 == LengthSequence-1);
		  while((l1 < (LengthSequence - OVERLAP)) || (l1 == 0))
			{
			  /** B.1. Measure G+C content in the current fragment: l1,l2 **/
			  GCScan(Sequence, GCInfo, l1, l2); 
			  GCScan(RSequence, GCInfo_r, LengthSequence-1 - l2,LengthSequence-1 - l1); 

			  /* G+C range: from 0 (l1) to l2 -l1 (l2) */
			  inigc = l1 -l1;
			  endgc = l2 -l1;
			  percentGC = ComputeGC(GCInfo,inigc,endgc); 
			  sprintf(mess,"G+C content in [%ld-%ld] is %f",l1, l2, percentGC);
			  printMess(mess); 
			  
			  /* Choose the isochore to predict sites according to the GC level */
			  currentIsochore = SelectIsochore(percentGC,isochores);
			  gp = isochores[currentIsochore];
			  sprintf(mess,"Selecting isochore %d", currentIsochore+COFFSET);
			  printMess(mess);

			  /* B.2. Prediction of sites and exons construction/filtering */
			  if (FWD) 
				{
				  /* Forward strand predictions */ 
				  sprintf(mess,"Running FWD  %s: %ld - %ld", Locus,l1,l2);
				  printMess(mess);
				  manager(Sequence, LengthSequence, 
						  allSites, allExons, l1, l2,
						  FORWARD, sr, gp,
						  isochores,nIsochores,
						  GCInfo);
				}      
			  if (RVS) 
				{
				  /* Reverse strand predictions */ 
				  sprintf(mess,"Running Reverse  %s: %ld - %ld(%ld - %ld)", 
						  Locus, LengthSequence-1 -l2, 
						  LengthSequence-1 -l1,l1,l2);         
				  printMess(mess);
				  manager(RSequence, LengthSequence,
						  allSites_r, allExons_r, 
						  LengthSequence-1 - l2, 
						  LengthSequence-1 - l1,
						  REVERSE, sr, gp,
						  isochores,nIsochores,
						  GCInfo_r);
				  				  
				  /* normalised positions: according to forward sense reading */
				  RecomputePositions(allSites_r, LengthSequence);
				  
				  /* exchange acc and donor sites to preserve Acc < Don */
				  SwitchPositions(allExons_r);
				}

			  /* B.3. Sort all of exons by left (minor) position */
			  if (EVD)
				{
				  /* Searching evidence exons in this split */
				  printMess("Searching annotations to be used in this split");
				  SearchEvidenceExons(evidence, (lastSplit)?l2:l2-OVERLAP);
				  
				  /* Unused annotations: out of range (info) */
				  if (lastSplit)
					{
					  sprintf(mess,"Forgotten last %ld evidences (out of range)",
							  evidence->nvExons - evidence->i2vExons);
					  printMess(mess);
					}
				}
			 
			  nExons = allExons->nExons + allExons_r->nExons;
			  if (EVD)
				{
				  nExons = nExons + evidence->ivExons;
				}

			  sprintf(mess,"Sorting %ld exons", nExons);  
			  printMess(mess);
			  
			  /* Merge predicted exons with some evidence exons */
			  SortExons(allExons, allExons_r, evidence, exons, l1, l2);
			  
			  /* Next block of annotations to be processed */
			  if (EVD)
				SwitchCounters(evidence);

			  /* B.4. Printing current fragment predictions (sites and exons) */
			  Output(allSites, allSites_r, allExons, allExons_r, 
					 exons, nExons, Locus, l1, l2, Sequence, gp, dAA); 

			  /* recompute stats about splice sites and exons */
			  updateTotals(m,allSites,allSites_r,allExons,allExons_r);
			  
			  /* B.5. Calling to genamic for assembling the best gene */ 
			  if (GENAMIC && nExons)
				{
				  genamic(exons, nExons, genes, gp);
				  
				  if (LengthSequence > LENGTHSi)
                    {
                      /* clean hash table of exons */
                      cleanDumpHash(dumpster->h);
                    } 

				  /* B.6. Backup operations of genes for the next split */
				  if (!lastSplit)
					{
					  /* backup of unused genes */
					  printMess("Back-up of d-genes");
					  BackupArrayD(genes, l2 - OVERLAP, gp, dumpster);
					  
					  /* back-up best partial genes */
					  printMess("Back-up of best partial genes\n");
					  BackupGenes(genes, gp->nclass, dumpster);
					}
				}
			  /* Computing new boundaries: next fragment in current sequence */
			  l1 += LENGTHSi - OVERLAP;
			  l2 = MIN(l1 + LENGTHSi -1, LengthSequence-1);
			  lastSplit = (l2 == LengthSequence-1);
			} /* processing next fragment */
		  
		  /* A.6. Full sequence processed: displaying best predicted gene */
		  if (GENAMIC)
			{	      
			  /* Printing gene predictions */
			  OutputGene(genes,
						 (EVD)? m->totalExons + evidence->nvExons : m->totalExons, 
						 Locus, Sequence, gp, dAA);

			  /* Reset best genes data structures for next input sequence */
			  printMess("Cleanning gene structures and dumpster");
			  cleanGenes(genes,gp->nclass,dumpster);
			}
		  
		  /* showing global stats about last sequence predicted */
		  OutputStats(Locus);

		  /* Reset evidence temporary counters */
		  if (EVD)
			resetEvidenceCounters(evidence);

		  cleanAcc(m);
		  strcpy(Locus,nextLocus);
		} /* endwhile(reading): next sequence to be processed... */
	} /*endifgeneid*/
  else
    {
      /* B. Only assembling genes from input exons */
      /* It might be necessary to increase NUMSITES and NUMEXONS */
	  
      /* B.0. Reading DNA sequence to make the translations */
      /* open the Sequence File */
      if ((seqfile = fopen(SequenceFile, "rb"))==NULL) 
        printError("The Sequence file can not be open for read");
      
      printMess("Reading DNA sequence");
      reading = IniReadSequence(seqfile,Locus);
      if (reading != EOF)
		{
		  reading = ReadSequence(seqfile, Sequence, nextLocus);
		  LengthSequence = FetchSequence(Sequence, RSequence);
		}
	  
      /* Header Output */
      OutputHeader(Locus, LengthSequence);
      
      /* B.1. Reading exons in GFF format */
      printMess("Reading exonsGFF from file");
      nExons = ReadExonsGFF(ExonsFile, evidence, isochores[0]->D);     
      sprintf(mess,"%ld exons GFF stored\n", nExons); 
      printMess(mess); 
      
      if (nExons) 
		{ 
		  /* B.2. Calling to genamic for assembling the best gene */
		  genamic(evidence->vExons, nExons, genes, isochores[0]);  
		}
      
      /* B.3. Printing gene predictions */
      OutputGene(genes, nExons, Locus, Sequence, isochores[0], dAA);
    } /* end only gene assembling from exons file */
  

  /* 4. The End */
  OutputTime();  
  
  exit(0);
  return(0);
}
