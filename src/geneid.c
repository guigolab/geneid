/*************************************************************************
*                                                                        *
*   Module: geneid                                                       *
*                                                                        *
*   geneid main program                                                  *
*                                                                        *
*   This file is part of the geneid distribution                         *
*                                                                        *
*     Copyright (C) 2006 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
*                          Tyler   ALIOTO                                * 
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


#include "geneid.h"
#include "bigbed.h"
#include "bigwig.h"
#ifdef WITH_HTSLIB
#include "bamcov.h"
#endif
/* #include <mcheck.h> */

/* geneid setup flags */
int
  /* sites to print */
  SFP=0, SDP=0, SAP=0, STP=0,
  /* exons to print */
  EFP=0, EIP=0, ETP=0, EXP=0, ESP=0, EOP = 0,
  /* introns to print */
  PRINTINT = 0,
  /* Partial or full prediction engine */
  GENAMIC = 1, GENEID = 1, 
  /* Only forward or reverse prediction engine */
  FWD=1, RVS=1,
  /* switch ORF prediction on */
  scanORF = 0,
  /* Input annotations or homology to protein information/reads to UTR prediction */
  EVD = 0, SRP = 0, UTR=0,
  /* Annotation-scoring mode (-J): score the forced evidence splice sites under
     the enabled profiles and classify each intron (report-only, no DP change) */
  SCOREANNOT = 0,
  /* -y library strandedness for BAM -S coverage: 0=unstranded (default),
     1=RF/dUTP, 2=FR (see BAMLIB_* in bamcov.h). Only used in a WITH_HTSLIB build. */
  BAMSTRAND = 0,
  /* Output formats */
  GFF = 0, GFF3 = 0, X10 = 0, XML = 0, cDNA = 0, PSEQ = 0, tDNA = 0,
  /* Verbose flag (memory/processing information) */
  BEG=0, VRB=0,
  /* Score for regions not-supported by protein homology */
  NO_SCORE, 
  /* Force single prediction: 1 gene */
  SGE=0,
  /* Detection of PolyPTracts in Acceptors */
  PPT=0,
  /* Detection of BranchPoints in Acceptors */
  BP=0,
  /* Detection of recursive splice sites */
  RSS=0,
  /* Detection of U12 introns */
  U12=0,
  /* Detection of U12gtag sites (acceptor uses BranchPoint)*/
  U12GTAG=0,
  /* Detection of U12atac sites (acceptor uses BranchPoint)*/
  U12ATAC=0,
  /* Detection of U2gcag sites */
  U2GCAG=0,
  /* Detection of U2gta donor sites */
  U2GTA=0,
  /* Detection of U2gtg donor sites */
  U2GTG=0,
  /* Detection of U2gty donor sites */
  U2GTY=0,
  /* Detection of PolyA Signal */
  PAS=0,
/* Length of flank around exons to subtract background RNA-seq signal */
  BKGD_SUBTRACT_FLANK_LENGTH = 0;


short
  /* Splice classes: the number of compatible splice site combinations used in genamic for joining exons */
  SPLICECLASSES = 1;
 
long
  /* User defined lower limit */
  LOW=0,
  /* User defined upper limit */
  HI=0;

float
  /* Millions of reads mapped */
MRM=15.0;

/* Optional Predicted Gene Prefix */
char  GenePrefix[MAXSTRING]="";
  

/* Increase/decrease exon weight value (exon score) */
float EW = NOVALUE;
float U12EW = 0; 
float EvidenceEW = 0; 
float EvidenceFactor = 1;
float U12_SPLICE_SCORE_THRESH = -1000;
float U12_EXON_SCORE_THRESH = -1000;

/* Weight applied to the U2 branch-point score when it is added to the acceptor
   score (see BuildAcceptors.c). Default 1 preserves the historical behaviour; a
   param file may set it to 0 so the branch is scored and reported (bp_score /
   bp_pos in the GFF) without contributing to the splice-site score. */
float BRANCH_SCORE_WEIGHT = 1;

/* Soft intron-length model (log-normal over ln(intron length)) emitted by
   geneid-train, plus the weight (lambda) on the smooth length-dependent penalty
   it drives on intron-spanning joins in genamic -- a soft replacement for the
   hard gene-model max distance. lambda default 0 leaves the penalty OFF, so the
   mu/sigma below are inert unless a param sets Intron_length_score_weight > 0. */
float INTRON_LENGTH_MU = 0;
float INTRON_LENGTH_SIGMA = 0;
float INTRON_LENGTH_WEIGHT = 0;

/* Detection of recursive splice sites */
float RSSMARKOVSCORE = 0;
float RSSDON = RDT;
float RSSACC = RAT;				  

/* Generic maximum values: sites, exons and backup elements */
long NUMSITES,NUMEXONS,MAXBACKUPSITES,MAXBACKUPEXONS,NUMU12SITES,NUMU12EXONS,NUMU12U12EXONS;

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
  long upperlimit;
  long lowerlimit;
  /* Forward semse data structures */
  packSites* allSites;
  packExons* allExons;
  
  /* Reverse sense data structures */
  packSites* allSites_r;
  packExons* allExons_r;
  
  /* Growable scratch buffers for sorting sites (see packSortSites) */
  packSortSites sortSites;

  /* Table to sort predicted exons by acceptor (growable; see SortExons) */
  exonGFF* exons;
  long exonscap;
  long nExons;
  
  /* External information: reannotation */
  packExternalInformation* external;
  packEvidence* evidence;
  packHSP* hsp;
  
  /* Best partial predicted genes */
  packGenes* genes;   
  
  /* Dumpster for backup operations between splits */
  packDump* dumpster;
  
  /* Amino acid dictionary (genetic code) */
  dict* dAA;
  
  /* geneid prediction parameters: data structures */
  gparam* gp = NULL;
  gparam** isochores;
 

    
  /* Input Filenames */
  char  SequenceFile[FILENAMELENGTH],
    ExonsFile[FILENAMELENGTH],
    HSPFile[FILENAMELENGTH],
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
  int reading;
  int lastSplit;
  BigBed* evBB = NULL;    /* non-NULL => -R evidence is a bigBed, queried per split */
  BamCov* evBam = NULL;   /* non-NULL => -R evidence is a BAM (introns), per split */
  long bbOwnedLo = 0;     /* upper acceptor bound owned by the previous fragment */
  char mess[MAXSTRING];

  
  /* Start memory trace -- for debugging memory leaks */
  /* mtrace(); */

  /** 0. Starting and reading options, parameters and sequence... **/
  nExons = 0;
  evidence = NULL;
  hsp = NULL;
  
  /* 0.a. Previous checkpoint about length in splits and overlapping */
  if (LENGTHSi <= OVERLAP)
    printError("LENGTHSi must be greater than OVERLAP parameter (geneid.h)");
  
  /* 0.b. Initializing stats and time counters */
  m = (account*)InitAcc();  
  
  /* 0.c. Read setup options */
  readargv(argc,argv,ParamFile,SequenceFile,ExonsFile,HSPFile,GenePrefix);
  printRes("\n\n\t\t\t** Running " GENEID_RELEASE " geneid@crg.es **\n\n");

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

  /* Estimation of memory required to execute geneid */
  if (BEG)
    beggar(LengthSequence);
  /** 1. Allocating main geneid data structures **/
  printMess("Request Memory to Operating System\n");
  
  /* 1.a. Mandatory geneid data structures */
  printMess("Request Memory Sequence\n");
  Sequence      = (char*)         RequestMemorySequence(LengthSequence);
  RSequence     = (char*)         RequestMemorySequence(LengthSequence);
  printMess("Request Memory Sites\n");
  allSites      = (packSites*)    RequestMemorySites();
  allSites_r    = (packSites*)    RequestMemorySites();
  printMess("Request Memory Exons\n");
  allExons      = (packExons*)    RequestMemoryExons();
  allExons_r    = (packExons*)    RequestMemoryExons();
  printMess("Request Memory Sort Exons\n");
  
  exons         = (exonGFF*)      RequestMemorySortExons();
  exonscap      = INITSORT;
  printMess("Request Memory Sort Sites\n");
  /* Scratch buffers for SortSites; each grows on demand (see GrowSiteArray) */
  sortSites.donorsites    = (site*) RequestMemorySortSites();
  sortSites.donorsitescap = INITSITESORT;
  sortSites.acceptorsites    = (site*) RequestMemorySortSites();
  sortSites.acceptorsitescap = INITSITESORT;
  sortSites.tssites = NULL;
  sortSites.tssitescap = 0;
  sortSites.tesites = NULL;
  sortSites.tesitescap = 0;
  if (UTR){
    sortSites.tssites    = (site*) RequestMemorySortSites();
    sortSites.tssitescap = INITSITESORT;
    sortSites.tesites    = (site*) RequestMemorySortSites();
    sortSites.tesitescap = INITSITESORT;
  }
  printMess("Request Memory Isochores, etc.\n");
  isochores     = (gparam**)      RequestMemoryIsochoresParams(); 
  GCInfo        = (packGC*)       RequestMemoryGC();
  GCInfo_r      = (packGC*)       RequestMemoryGC();
  dAA           = (dict*)         RequestMemoryAaDictionary();
  printMess("Request Memory External\n");
  external      = (packExternalInformation*) RequestMemoryExternalInformation();
  

  /* 1.b. Backup information might be necessary between splits */
  if (LengthSequence > LENGTHSi){
    printMess("Request Memory Dumpster\n");
    dumpster   = (packDump*)     RequestMemoryDumpster();
  }else{
    dumpster = NULL;
  }
  /** 2. Reading statistical model parameters file **/
  printMess("Reading parameters..."); 
  nIsochores = readparam(ParamFile, isochores);
  if (U12){
    if ((!U12GTAG)&&(!U12ATAC)){
      U12GTAG = 0;
      U12ATAC = 0;
    }
  }else if (SCOREANNOT && !GENEID){
    /* -J -O annotation scoring: keep the param's U12 profiles so introns can be
       typed U12. Safe because -J scoring is report-only and runs AFTER genamic
       (no ab-initio U12 building happens in the assembly-only path), and -U --
       the usual U12 switch -- is disallowed alongside -O. */
  }else{
    U12GTAG = 0;
    U12ATAC = 0;
  }
  /** 1. Allocating genes data structure (after the number of splice classes has been determined **/
  genes = (packGenes*) RequestMemoryGenes();

  /** 3. Starting processing. This per-locus loop is SHARED between full
     prediction and assemble-only (-O): when GENEID is off, the ab-initio
     prediction (manager) is skipped and the exons come from the -O file via the
     evidence path, while the multi-locus loop, isochore selection, gene
     assembly and output are all inherited. **/
    {
      /* A. Predicting signals, exons and genes in DNA sequences */
      /* A.1. Reading external information I: annotations */
      if (EVD)
	{
	  /* bbOpen/bamOpen validate their file magic and return NULL for a text
	     GFF (whose first bytes are ASCII), so they double as the format sniff:
	     bigBed -> per-split records; BAM -> per-split spliced-read introns
	     (WITH_HTSLIB build); otherwise the text annotation path. */
	  evBB = bbOpen(ExonsFile);
	  if (evBB)
	    printMess("Reading evidence from bigBed (per-split range queries)...");
#ifdef WITH_HTSLIB
	  else if ((evBam = bamOpen(ExonsFile)) != NULL)
	    printMess("Reading introns from BAM junctions (per-split range queries)...");
#endif
	  else
	    {
	      printMess("Reading evidence (annotations)...");
	      external->nvExons =
		ReadExonsGFF(ExonsFile, external, isochores[0]->D);
	      sprintf(mess,"%ld annotations acquired from file\n",
		      external->nvExons);
	      printMess(mess);
	    }
	}
	
      /* A.2. Reading external information II: homology / RNA-seq coverage.
	 -S plus.bw,minus.bw = stranded bigWig coverage; -S cov.bw = unstranded
	 bigWig; -S reads.bam = indexed BAM coverage (WITH_HTSLIB build); otherwise
	 the text HSP path (ReadHSP). bwOpen/bamOpen validate their file magic, so
	 they double as the format sniff (a text HSP GFF matches neither). */
      if (SRP)
	{
	  char* comma = strchr(HSPFile, ',');
	  if (comma != NULL)
	    {
	      *comma = '\0';
	      external->bwPlus  = bwOpen(HSPFile);
	      external->bwMinus = bwOpen(comma + 1);
	      *comma = ',';
	      if (external->bwPlus == NULL || external->bwMinus == NULL)
		printError("-S with two comma-separated files expects stranded bigWigs (plus.bw,minus.bw)");
	    }
	  else
	    {
	      external->bwPlus = bwOpen(HSPFile);
	      external->bwMinus = external->bwPlus;   /* unstranded: same signal both strands */
#ifdef WITH_HTSLIB
	      if (external->bwPlus == NULL)
		external->bam = bamOpen(HSPFile);   /* not a bigWig: try an indexed BAM */
#endif
	    }

	  if (external->bam != NULL)
	    {
	      /* Coverage fills sr[] to score exons with or without -u; -u adds UTR
		 prediction (and the rpkm report) on top. */
	      printMess("Reading RNA-seq coverage from indexed BAM (per-split range queries)...");
	    }
	  else if (external->bwPlus != NULL)
	    {
	      printMess("Reading RNA-seq coverage from bigWig (per-split range queries)...");
	    }
	  else
	    {
	      printMess("Reading homology information...");
	      external->nHSPs = ReadHSP(HSPFile, external);
	      sprintf(mess,"%ld HSPs acquired from file",
		      external->nHSPs);
	      printMess(mess);
	    }
	}
	  
      if (EVD || SRP)
	{
	  sprintf(mess,"External information acquired from %ld sequences\n",
		  external->nSequences);
	  printMess(mess); 
	}

      /** A.3. Input DNA sequences (perhaps more than one) **/
      if ((seqfile = fopen(SequenceFile, "rb"))==NULL) 
        printError("The input sequence file can not be accessed");
	  
      /* reading the locusname of sequence (in Fasta format) */
      reading = IniReadSequence(seqfile,Locus);
	  
      while (reading != EOF)
	{		  
          printMess("Loading DNA sequence");
	  reading = ReadSequence(seqfile, Sequence, nextLocus);
		  		  
          /* A.3. Prepare sequence to work on */
          printMess("Processing DNA sequence");
          LengthSequence = FetchSequence(Sequence, RSequence);
	  OutputHeader(Locus, LengthSequence);

	  /* name of the current sequence, for the per-fragment bigWig query */
	  external->curLocus = Locus;

	  /* A.4. Prepare external information */
	  if (SRP && external->bwPlus == NULL)   /* text HSP path only */
	    {
	      printMess("Select homology information");
	      hsp = (packHSP*) SelectHSP(external, Locus, LengthSequence);
	      if (hsp == NULL)
		sprintf(mess,"No information has been provided for %s\n",
			Locus);
	      else
		sprintf(mess,"Using %ld HSPs in %s\n",
			hsp->nTotalSegments,Locus);

	      printMess(mess);
	    }

	  if (EVD)
	    {
	      if (evBB || evBam)
		/* bigBed / BAM: no preload -- evidence[0] is filled per fragment below */
		evidence = external->evidence[0];
	      else
		{
		  printMess("Select annotations");
		  evidence = (packEvidence*) SelectEvidence(external,Locus);
		  if (evidence == NULL)
		    sprintf(mess,"No information has been provided for %s\n",
			    Locus);
		  else
		    sprintf(mess,"Using %ld annotations in %s\n",
			    evidence->nvExons,Locus);
		  printMess(mess);
		}
	    }

	  /* A.5. Processing sequence into several fragments if required */
	  /* l1 is the left end and l2 is the right end in Sequence */
	  /* The arguments HI and LO are converted into lower and upper limit coordinates and l1 and l2 are adjusted */
	  /* -j/-k (LOW/HI) restrict prediction to a sub-window [lowerlimit,
	     upperlimit] of the sequence (defaulting to the whole thing); THAT
	     window is what actually gets split into LENGTHSi-sized, OVERLAP-
	     overlapping fragments [l1,l2] below -- each new fragment starts
	     OVERLAP bp before the previous one ended (l1 += LENGTHSi-OVERLAP),
	     so a splice signal near a fragment boundary is still seen with
	     full context by at least one fragment (see manager.c's boundary-
	     window comment for how that overlap is actually used). Everything
	     genamic assembles that might still be needed once this fragment's
	     arrays are reused gets copied into the dumpster first (see B.6
	     below and BackupGenes.c's file comment) -- that's the whole
	     reason the dumpster exists. */
	  upperlimit = LengthSequence-1;
	  lowerlimit = 0;

	  if ((HI > 0)&&(HI >= LOW)&&(HI < LengthSequence)){
	    upperlimit = HI - 1;
	  }else{
	    upperlimit = LengthSequence-1;
	  }
	  if ((LOW > 0)&&(LOW <= upperlimit)){
	    lowerlimit = LOW - 1;
	  }else{lowerlimit = 0;}
	  bbOwnedLo = lowerlimit;   /* first fragment owns acceptors above this */
	  l1 = lowerlimit;
	  l2 = MIN(l1 + LENGTHSi-1,LengthSequence-1);
	  l2 = MIN(l2,upperlimit);
	  /* Check to see if we are on last split */
	  lastSplit = (l2 == upperlimit);
	  sprintf(mess,"Running on range %ld to %ld\n",
		  lowerlimit ,upperlimit);
	  printMess(mess);
	  /* Runs at least once (l1==lowerlimit on the first pass) even if the
	     whole window is short enough to need only one fragment; otherwise
	     continues until l1 has advanced within OVERLAP of upperlimit. */
	  while((l1 < (upperlimit + 1 - OVERLAP)) || (l1 == 0)|| (l1 == lowerlimit))
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

	      /* B.2. Prediction of sites and exons construction/filtering.
		 Skipped entirely when GENEID is off (-O): the exons then come only
		 from the evidence/-O file, merged below by SortExons. */
	      if (GENEID && FWD)
		{
		  /* Forward strand predictions */
		  sprintf(mess,"Running FWD  %s: %ld - %ld", Locus,l1,l2);
		  printMess(mess);
		  manager(Sequence, LengthSequence, 
			  allSites, allExons, 
			  l1, l2, lowerlimit, upperlimit,
			  FORWARD, 
			  external, hsp, gp,
			  isochores,nIsochores,
			  GCInfo,&sortSites);
		}      
	      if (GENEID && RVS)
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
			  LengthSequence-1 - upperlimit, LengthSequence-1 - lowerlimit,
			  REVERSE, 
			  external, hsp, gp,
			  isochores,nIsochores,
			  GCInfo_r,&sortSites);
				  				  
		  /* normalised positions: according to forward sense reading */
		  RecomputePositions(allSites_r, LengthSequence);
				  
		  /* exchange acc and donor sites to preserve Acc < Don */
		  SwitchPositions(allExons_r);
		}

	      /* B.3. Sort all of exons by left (minor) position */
	      if (EVD && evidence != NULL)
		{
		  /* Searching evidence exons in this fragment */
		  printMess("Searching annotations to be used in this fragment");
		  if (evBB || evBam)
		    {
		      /* Per-split bigBed / BAM query. ownedHi mirrors SearchEvidenceExons'
			 (l2 or l2-OVERLAP) bound in 1-based acceptor units, so each
			 record is committed in exactly one fragment. */
		      long ownedHi = ((lastSplit)?l2:l2-OVERLAP) + 1;
		      if (evBB)
			ReadExonsBigBed(evBB, external, isochores[0]->D,
					Locus, l1, l2, bbOwnedLo, ownedHi);
#ifdef WITH_HTSLIB
		      else
			ReadIntronsBam(evBam, external, isochores[0]->D,
				       Locus, l1, l2, bbOwnedLo, ownedHi,
				       Sequence, LengthSequence);
#endif
		      bbOwnedLo = ownedHi;
		    }
		  else
		    {
		      SearchEvidenceExons(external,
					  evidence,
					  (lastSplit)?l2:l2-OVERLAP);

		      /* Unused annotations: out of range (info) */
		      if (lastSplit)
			{
			  sprintf(mess,"Leaving out last %ld evidences (out of range)",
				  evidence->nvExons - external->i2vExons);
			  printMess(mess);
			}
		    }
		}
			 
	      nExons = allExons->nExons + allExons_r->nExons;
	      if (EVD && evidence != NULL)
		nExons = nExons + external->ivExons;

	      /* BEGIN artificial exon: + and - */
	      if (l1 == lowerlimit){
		nExons = nExons + 2;
	      }
	      /* END artitificial exon: + and - */
	      if (l2 == upperlimit){
		nExons = nExons + 2;
	      }
/* 	      sprintf(mess,"l1: %ld   ll:%ld   l2: %ld   ul: %ld\n", l1,lowerlimit,l2,upperlimit); */
/* 	      printMess(mess); */
/* 	      /\* B.4. Printing current fragment predictions (sites and exons) *\/ */
/* 	      Output(allSites, allSites_r, allExons, allExons_r,  */
/* 		     exons, nExons, Locus, l1, l2, lowerlimit, Sequence, gp, dAA, GenePrefix);  */

	      sprintf(mess,"Sorting %ld exons\n", nExons);  
	      printMess(mess);
			  
	      /* Merge predicted exons with some evidence exons */
	      SortExons(allExons, allExons_r,
			external,
			evidence,
			&exons, &exonscap,
			l1, l2,
			lowerlimit,
			upperlimit);
	      sprintf(mess,"Finished sorting %ld exons\n", nExons);  
	      printMess(mess);
			  
	      /* Next block of annotations to be processed (GFF cursor only; the
		 bigBed/BAM paths re-set the window per fragment). */
	      if (EVD && evidence != NULL && !evBB && !evBam)
		SwitchCounters(external);

	      /* B.4. Printing current fragment predictions (sites and exons) */
	      Output(allSites, allSites_r, allExons, allExons_r, 
		     exons, nExons, Locus, l1, l2, lowerlimit, Sequence, gp, dAA, GenePrefix); 

	      /* recompute stats about splice sites and exons */
	      updateTotals(m,allSites,allSites_r,allExons,allExons_r);
			  
	      /* B.5. Calling to genamic for assembling the best gene */ 
	      if (GENAMIC && nExons)
		{

		  genamic(exons, nExons, genes, gp);
				  
		  if (upperlimit - lowerlimit + 1 > LENGTHSi)/*  if (LengthSequence > LENGTHSi) */
		    {
		      /* clean hash table of exons */
		      cleanDumpHash(dumpster->h);
		    } 

		  /* B.6. Backup operations of genes for the next split */
		  if (!lastSplit)
		    {
		      /* backup of unused genes */
		      printMess("Back-up of d-genes");
		      /* l2-OVERLAP: the next fragment starts at l1+LENGTHSi-OVERLAP
			 (<= l2), so exons ending before l2-OVERLAP are behind where
			 the next fragment will (re)scan from and can be trimmed --
			 see BackupArrayD's own comment in BackupGenes.c. */
		      BackupArrayD(genes, l2 - OVERLAP, gp, dumpster);

		      /* back-up best partial genes */
		      printMess("Back-up of best partial genes\n");
		      BackupGenes(genes, gp->nclass, dumpster);
		    }
		}
	      /* Computing new boundaries: next fragment in current sequence */
	      l1 += LENGTHSi - OVERLAP;
	      l2 = MIN(l1 + LENGTHSi -1, upperlimit);
	      lastSplit = (l2 == upperlimit);
	    } /* processing next fragment */
		  
	  /* Annotation-scoring mode (-J): score + classify the evidence splice sites
	     of the assembled best gene (genes->GOptim, the same chain OutputGene is
	     about to print) for reporting. Runs AFTER all fragments are assembled, so
	     it is strictly report-only (never affects the DP). Scoring the PRINTED
	     chain -- not the loaded evidence -- is what keeps it correct on multi-split
	     sequences, where the printed exons are dumpster deep-copies. Single
	     isochore (isochores[0]). */
	  if (SCOREANNOT && EVD && evidence != NULL)
	    ScoreEvidenceSites(genes->GOptim, Sequence, RSequence, LengthSequence, isochores[0]);

	  /* A.6. Full sequence processed: displaying best predicted gene */
	  if (GENAMIC)
	    {
	      /* Printing gene predictions */
	      OutputGene(genes,
			 (EVD && evidence != NULL)? 
			 m->totalExons + evidence->nvExons : 
			 m->totalExons, 
			 Locus, Sequence, gp, dAA, GenePrefix);

	      /* Reset best genes data structures for next input sequence */
	      printMess("Cleaning gene structures and dumpster");
	      cleanGenes(genes,gp->nclass,dumpster);
	    }
		  
	  /* showing global stats about last sequence predicted */
	  OutputStats(Locus);

	  /* Reset evidence temporary counters */
	  if (EVD && evidence != NULL)
	    resetEvidenceCounters(external);

	  cleanAcc(m);
	  strcpy(Locus,nextLocus);
	} /* endwhile(reading): next sequence to be processed... */
    } /* end shared per-locus processing (prediction and -O assemble-only) */
  

  /* Close bigWig coverage readers (if any). When unstranded, bwPlus==bwMinus. */
  if (external->bwMinus != NULL && external->bwMinus != external->bwPlus)
    bwClose(external->bwMinus);
  if (external->bwPlus != NULL)
    bwClose(external->bwPlus);
#ifdef WITH_HTSLIB
  if (external->bam != NULL)
    bamClose(external->bam);
  if (evBam != NULL)
    bamClose(evBam);
#endif

  /* 4. The End */
  OutputTime();

  exit(0);
  return(0);
}
