/*************************************************************************
*                                                                        *
*   Module: readargv                                                     *
*                                                                        *
*   Read set up options and filenames from user input                    *
*                                                                        *
*   This file is part of the geneid distribution                         *
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


#include "geneid.h"

/* geneid.c external vars */
extern int  SFP,SDP,SAP,STP,
            EFP,EIP,ETP,EXP,ESP,EOP,
            U12, PRINTINT,RSS,
            VRB,
            FWD,RVS,
            GENEID, GENAMIC,
            GFF, GFF3, X10,
            EVD, SRP, BEG, UTR,
            scanORF, XML, cDNA, PSEQ, tDNA,
            SGE, SCOREANNOT;
extern float EW,EvidenceEW,MRM;
extern int MRMset;   /* set here when -N is given, so a BAM's index is not used to estimate MRM */
extern int EXPRLLR;        /* -L: enable Poisson per-base expression LLR coverage scoring */
extern float LLRK, LLRW;   /* -L fold-change k (>1); -Q weight/scale of the LLR term */
extern int LLRWset;        /* set by -Q: suppresses the auto-derived weight */
extern long LOW,HI;
extern int BAMSTRAND;   /* -y library strandedness for BAM -S coverage */

/* required by getopts */
extern char* optarg;
extern int optind;

char* USAGE="NAME\n\tgeneid - a program to annotate genomic sequences\nSYNOPSIS\n\tgeneid\t[-bdaefitnxszru]\n\t\t[-TDAZU]\n\t\t[-p gene_prefix]\n\t\t[-G] [-3] [-X] [-M] [-m]\n\t\t[-WCF] [-o] [-J]\n\t\t[-j lower_bound_coord]\n\t\t[-k upper_bound_coord]\n\t\t[-N numer_nt_mapped]\n\t\t[-L fold_change] [-Q llr_weight]\n\t\t[-O <gff_exons_file>]\n\t\t[-R <gff_annotation-file>]\n\t\t[-S <gff_homology_file>]\n\t\t[-Y reads.bam]\n\t\t[-y library_type]\n\t\t[-P <parameter_file>]\n\t\t[-E exonweight]\n\t\t[-V evidence_exonweight]\n\t\t[-Bv] [-h]\n\t\t<locus_seq_in_fasta_format>\nRELEASE\n\t" GENEID_RELEASE "\n";

void printHelp()
{
  printf("%s", USAGE);

  printf ("OPTIONS\n");
  
  printf("\t-b: Output Start codons\n");
  printf("\t-d: Output Donor splice sites\n");
  printf("\t-a: Output Acceptor splice sites\n");
  printf("\t-e: Output Stop codons\n");
  
  printf("\t-f: Output Initial exons\n");
  printf("\t-i: Output Internal exons\n");
  printf("\t-t: Output Terminal exons\n");
  printf("\t-n: Output introns\n");
  printf("\t-s: Output Single genes\n");
  printf("\t-x: Output all predicted exons\n");
  printf("\t-z: Output Open Reading Frames\n\n");
  
  printf("\t-T: Output genomic sequence of exons in predicted genes\n");
  printf("\t-D: Output genomic sequence of CDS in predicted genes\n");
  printf("\t-A: Output amino acid sequence derived from predicted CDS\n\n");
  printf("\t-p: Prefix this value to the names of predicted genes, peptides and CDS\n\n");

  printf("\t-G: Use GFF format to print predictions\n");
  printf("\t-3: Use GFF3 format to print predictions\n");
  printf("\t-X: Use extended-format to print gene predictions\n");
  printf("\t-M: Use XML format to print gene predictions\n");
  printf("\t-m: Show DTD for XML-format output \n\n");

  printf("\t-j  <coord>: Begin prediction at this coordinate\n");
  printf("\t-k  <coord>: End prediction at this coordinate\n");  
  printf("\t-N  <num_reads>: Millions of reads mapped to genome (rpkm report; for a\n"
	 "\t     BAM input this is estimated from the index when -N is omitted)\n");
  printf("\t-L  <k>: enable Poisson expression LLR coverage scoring; k>1 is the\n"
	 "\t     enrichment of an expressed exon over the genome-wide mean coverage\n"
	 "\t     (tens, not units). Default off (legacy log(cov+1)/raw term).\n"
	 "\t     Recommended start: -L 50 -Q 0.0007 (human RNA-seq, bam+u)\n");
  printf("\t-Q  <w>: weight/scale of the -L LLR term (default 0.0007). Transfers\n"
	 "\t     across library depths; re-tune per param file\n");
  printf("\t-W: Only Forward sense prediction (Watson)\n");
  printf("\t-C: Only Reverse sense prediction (Crick)\n");
  printf("\t-U: Allow U12 introns (Requires appropriate U12 parameters to be set in the parameter file)\n");
  printf("\t-r: Use recursive splicing\n");
  printf("\t-F: Force the prediction of one gene structure\n");
  printf("\t-o: Only running exon prediction (disable gene prediction)\n");
  printf("\t-O  <exons_filename>: Only running gene prediction (not exon prediction)\n");
  printf("\t-Z: Activate Open Reading Frames searching\n\n");
  
  printf("\t-R  <exons_filename>: Provide annotations to improve predictions (GFF, or a bigBed of the same records)\n");
  printf("\t    In a WITH_HTSLIB build, -R may instead be an indexed BAM: spliced-read (CIGAR N) junctions become Intron evidence\n");
  printf("\t    (junction strand from the XS tag, else minimap2 ts, else the GT-AG/CT-AC splice motif)\n");
  printf("\t-J: Annotation-scoring mode: score the splice sites of forced annotation\n");
  printf("\t    exons (from -O or -R) under the parameter's profiles and classify each\n");
  printf("\t    intron as U2 or U12 (report-only; does not change the assembly). Best\n");
  printf("\t    used as -J -O <annotation> to score/type a provided gene structure\n");
  printf("\t-S  <HSP_filename>: Using information from protein sequence alignments to improve predictions\n");
  printf("\t    RNA-seq coverage may instead be given as bigWig(s): -S plus.bw,minus.bw (stranded) or -S cov.bw (unstranded)\n");
  printf("\t    or, in a WITH_HTSLIB build, as an indexed BAM: -S reads.bam. Coverage scores exons with or without -u (-u adds UTR prediction)\n");
  printf("\t-Y  <reads.bam>: (WITH_HTSLIB) one indexed BAM as BOTH intron evidence (-R junctions) and RNA-seq coverage (-S)\n");
  printf("\t-y  <rf|fr|none>: BAM -S library strandedness -- rf=dUTP/reverse, fr=forward, none=unstranded (default)\n\n");
  printf("\t-u: Turn on UTR prediction. Only valid with -S option: HSP/EST/short read ends are used to determine UTR ends\n");
  
  printf("\t-E: Add this value to the exon weight parameter (see parameter file)\n");
  printf("\t-V: Add this value to the score of evidence exons \n");
  printf("\t-P  <parameter_file>: Use other than default parameter file (human)\n\n");
  
  printf("\t-B: Display memory required to execute geneid given a sequence\n");
  printf("\t-v: Verbose. Display info messages\n");
  printf("\t-h: Show this help\n");

  printf ("AUTHORS\n");
  printf("\t%s has been developed by Enrique Blanco, Tyler Alioto and Roderic Guigo.\n\tParameter files have been created by Genis Parra and Tyler Alioto. Any bug or suggestion\n\tcan be reported to geneid@crg.es\n",GENEID_RELEASE);

  printf("\n\n\n");
}

void printDTD()
{
  printf("<?xml version=\"1.0\" ?>");
  printf("<!-- DTD for XML format in geneid output -->");

   printf("<!-- Element declarations -->");
   printf("<!ELEMENT prediction (gene*)>");
   printf("<!ELEMENT gene ((exon+),cDNA,protein)>");
   printf("<!ELEMENT exon (site,site)>");
   printf("<!ELEMENT cDNA (#PCDATA)>");
   printf("<!ELEMENT protein (#PCDATA)>");

   printf("<!-- Attribute declarations -->");
   printf("<!ATTLIST prediction");
        printf("\tlocus    CDATA   #REQUIRED");
        printf("\tlength   CDATA   #IMPLIED");
        printf("\tsource   CDATA   #IMPLIED");
        printf("\tdate     CDATA   #IMPLIED");
        printf("\tgenes    CDATA   #REQUIRED");
	    printf("\tscore    CDATA   #REQUIRED>");

   printf("<!ATTLIST gene");
        printf("\tidGene   ID      #REQUIRED");
        printf("\tstrand   (fwd|rvs)   #IMPLIED");
        printf("\tnExons   CDATA   #IMPLIED");
        printf("\tscore    CDATA   #REQUIRED>");

   printf("<!ATTLIST exon");  
        printf("\tidExon   ID      #REQUIRED");
        printf("\ttype     (First | Internal | Terminal | Single) #REQUIRED");
	printf("frame    (0|1|2) #REQUIRED");
        printf("\tscore    CDATA   #REQUIRED>");

   printf("<!ATTLIST site");
        printf("\tidSite   ID      #REQUIRED");
        printf("\ttype     (Acceptor | Donor | Start | Stop) #REQUIRED");
	printf("position CDATA   #REQUIRED");
        printf("\tscore    CDATA   #REQUIRED>\n\n");
}

/* Parses the command line via getopt into the external flags declared
 * above (shared with geneid.c/manager.c), then validates the combination.
 * The pipeline has two stages that can each be switched off: GENEID (the
 * ab initio exon-prediction stage, turned off by -O <exons_file>, which
 * supplies exons directly instead) and GENAMIC (the gene-assembly/DP
 * stage, turned off by -o, exon-prediction-only mode -- see Output.c's
 * S*P/E*P flags for what gets printed instead). geneidOpts/genamicOpts
 * tally how many options were given that only make sense when the
 * corresponding stage actually runs, so the checks below can reject e.g.
 * -O combined with an exon-only-stage option. geneidOpts counts only the
 * PREDICTION-ONLY options -- ones with no effect once ab initio prediction is
 * off (-R/-S evidence & homology into the scorer, -Z ORF scan, -r recursive
 * splice, and the -a/-b/-d/-e/-f/-i/-s/-t/-x/-z predicted-feature print flags).
 * Options that also make sense on the assemble-only path (-u UTR, -U U12 typing,
 * -N reads-mapped, -E exon weight, -F single-gene) are deliberately NOT counted,
 * so they are allowed alongside -O (assemble + score a provided annotation).
 * printOptions tallies the
 * single-feature debug print flags (-b/-d/-a/-e/-f/-i/-t/-s/-x/-z, i.e.
 * SFP/SDP/SAP/STP/EFP/EIP/ETP/ESP/EXP/EOP), which are mutually exclusive
 * with XML output (-M) below. */
void readargv (int argc,char* argv[],
			   char* ParamFile, char* SequenceFile,
	       char* ExonsFile, char* HSPFile, char* GenePrefix)
{
  int c;
  int error=0;
  int geneidOpts = 0;
  int genamicOpts = 0;
  int printOptions =0;
  int NOpt =0;
  char mess[MAXSTRING];
  char *dummy1;
  char *dummy2;
  char *dummy3;
  /* Reading setup options */
  while ((c = getopt(argc,argv,"oO:bdaefitnsrxj:k:N:p:UDATzZXmMG3BvE:V:R:S:WCFP:huJy:Y:L:Q:")) != -1)
    switch(c)
      {
      case 'B': BEG++; 
		break;
	  case 'C': FWD--;
		/* geneidOpts++; */
		break;
	  case 'T': tDNA++;
		genamicOpts++;
		break;
	  case 'D': cDNA++;
		genamicOpts++;
		break;
	  case 'A': PSEQ++;
		genamicOpts++;
		break;
	  case 'p': strcpy (GenePrefix,optarg);
		genamicOpts++;
		break;
	  case 'E': EW = atof(optarg);
		/* assembly-compatible: exon weight, allowed under -O */
		break;
	  case 'V': EvidenceEW = atof(optarg);
		genamicOpts++;
		break;
	  case 'F': SGE++;
		/* assembly-compatible: force single gene, allowed under -O */
		break;
	  case 'G': GFF++;
		break;
	  case '3': GFF3++;
	  	GFF++;
		break;
	  case 'M': XML++;
	        genamicOpts++;
		break;
	  case 'O': GENEID--;   /* assemble only: skip prediction, feed exons as evidence */
		EVD++;
		strcpy (ExonsFile,optarg);
		break;
	  case 'P': strcpy (ParamFile,optarg); 
		break;
          case 'R': EVD++;
		strcpy (ExonsFile,optarg);
		geneidOpts++;
		break;
          case 'S': SRP++;
		strcpy (HSPFile,optarg);
		geneidOpts++;
		break;
	  case 'Y':   /* one indexed BAM as BOTH intron evidence (-R) and RNA-seq coverage (-S) */
#ifndef WITH_HTSLIB
		printError("-Y (one BAM for introns + coverage) requires building geneid with WITH_HTSLIB=1");
#endif
		EVD++;
		strcpy (ExonsFile,optarg);
		SRP++;
		strcpy (HSPFile,optarg);
		geneidOpts++;
		break;
          case 'u': UTR++;
		/* assembly-compatible: assemble UTR exons from the annotation, allowed under -O */
		break;
	  case 'W': RVS--;
		/* geneidOpts++; */
		break; 
          case 'X': X10++;
		genamicOpts++;
		break;
          case 'Z': scanORF++;
		geneidOpts++;
		break;      
	  case 'a': SAP++;
		geneidOpts++;
		printOptions++;
		break;
          case 'b': SFP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'd': SDP++;
		geneidOpts++;
		printOptions++;
		break;
	  case 'e': STP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'f': EFP++;
		geneidOpts++;
		printOptions++;
		break;
	  case 'h': printHelp();
		exit(0);
		break;
      case 'i': EIP++;
		geneidOpts++;
		printOptions++;
		break;
	  case 'j': LOW = strtol(optarg,&dummy1,0);
		/* geneidOpts++; */
		break;
	  case 'k': HI = strtol(optarg,&dummy2,0);
		/* geneidOpts++; */
		break;
	  case 'N': MRM = strtof(optarg,&dummy3);
		/* assembly-compatible: reads-mapped (rpkm reporting), allowed under -O */
		MRMset = 1;   /* explicit value: suppress BAM-index auto-estimation */
		NOpt++;
		break;
	  case 'L': EXPRLLR = 1;   /* enable Poisson expression LLR coverage scoring */
		LLRK = strtof(optarg,&dummy3);
		if (LLRK <= 1.0)
		  printError("-L expects a fold-change k > 1 (expressed/background)");
		break;
	  case 'Q': LLRW = strtof(optarg,&dummy3);   /* LLR weight/scale */
		LLRWset = 1;   /* explicit: suppress the auto-derived weight */
		break;
	  case 'U': U12++;
		/* assembly-compatible: U12 intron typing, allowed under -O */
		break;
	  case 'r': RSS++;
		geneidOpts++;
		break;
	  case 'm': printDTD();
		exit(0);
		break;
	  case 'o': GENAMIC--;
		break;
      case 's': ESP++;
		geneidOpts++;
		printOptions++;
		break;
	  case 't': ETP++;
		printOptions++;
		geneidOpts++;
		break;
	  case 'n': PRINTINT++;
		printOptions++;
		/* geneidOpts++; */
		break;
	  case 'J': SCOREANNOT++;   /* score forced-annotation splice sites + type introns */
		break;
	  case 'y':   /* library strandedness for BAM -S coverage */
		if (!strcmp(optarg,"rf") || !strcmp(optarg,"RF"))
		  BAMSTRAND = 1;
		else if (!strcmp(optarg,"fr") || !strcmp(optarg,"FR"))
		  BAMSTRAND = 2;
		else if (!strcmp(optarg,"none") || !strcmp(optarg,"unstranded"))
		  BAMSTRAND = 0;
		else
		  printError("-y expects a library type: rf (dUTP), fr, or none");
		break;
      case 'v': VRB++;
		break;
	  case 'x': EXP++;
		geneidOpts++;
		printOptions++;
		break;
      case 'z': EOP++;
		geneidOpts++;
		printOptions++;
		break;
      }
  
  /* Setup Errors (a): Incompatible options selected */
  if (!GENEID && geneidOpts)
    printError("Incompatible options (with -O)");
  
  if (!GENAMIC && genamicOpts)
    printError("Incompatible options (with -o)");
  
  if (!GENAMIC && !GENEID)
    printError("Incompatible options (-o | -O)");
 
  if (XML && printOptions)
    printError("Incompatible options (-M | print gene features)"); 

  if (XML && (GFF || X10))
    printError("Incompatible options (XML and other output formats)");

  if (cDNA && GFF && !GFF3)
    printError("Incompatible options( -D | -G)");
  
  if (PSEQ && GFF && !GFF3)
    printError("Incompatible options( -A | -G)");

  /* UTR prediction from RNA-seq needs the -S coverage file, but only when ab
     initio prediction runs. Under -O (GENEID off) the UTR exons are supplied in
     the annotation and assembled directly, so -S is not required there. */
  if (UTR && !SRP && GENEID)
    printError("UTR option ( -u )selected without required SR file ( -S <HSP_filename>)");
  if (NOpt && !UTR)
    printError("N option requires UTR option ( -u )");
  if (error)
	{
	  sprintf(mess,"Wrong usage of options\n%s",USAGE);
	  printError(mess);
	}

  /* Setup Errors (b): Wrong number of filenames */
  /* The one non-option argument left after getopt is the input FASTA file,
     required unless -O already disabled GENEID (gene-assembly-only runs
     from a pre-supplied ExonsFile need no separate sequence input here). */
  /* Read the name of the input fasta file */
  if (optind < argc)
    {
      strcpy(SequenceFile,argv[optind]);
      optind++;
      if (optind < argc)
		{ 
          sprintf(mess,"Input contains more than one file but only one is required\n%s",USAGE);
		  printError(mess);
		}
    }
  else
    if (GENEID)
      {
        sprintf(mess,"One filename is required (DNA sequence, Fasta format)\n%s",USAGE);
		printError(mess);
      }
  
  /* Default parameter file selected if option -P not used */
  if (!strcmp(ParamFile,""))
    strcpy(ParamFile,PARAMETERFILE);
}
