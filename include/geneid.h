/*************************************************************************
*                                                                        *
*   Module: geneid.h                                                     *
*                                                                        *
*   Definitions, data structures types and imported headers              *
*                                                                        *
*   This file is part of the geneid 1.2 distribution                     *
*                                                                        *
*     Copyright (C) 2003 - Enrique BLANCO GARCIA                         *
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

/* $Id: geneid.h,v 1.54 2010/11/25 20:48:06 talioto Exp $ */

/* Required libraries */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

/*************************************************************************
A. DEFINITIONS
*************************************************************************/

/* The name of the game                     */
#define VERSION   "geneid_v1.4"  
#define SITES     "geneid_v1.4"  
#define EXONS     "geneid_v1.4"       
#define EVIDENCE  "evidence"           

/* -------------------------------------------------------------------------
 * These constants used to be a build-time "memory profile": geneid reserved
 * fixed-size arrays up front (NUMSITES/NUMEXONS x FSORT, the per-type build
 * arrays, the dumpster, ...), and the Makefile tuned them with -D to trade
 * memory for capacity. All of those arrays now GROW ON DEMAND, so there is no
 * profile to pick and nothing here is a correctness ceiling any more. What is
 * left are just derivation constants: LENGTHSi (the fragment window -- an
 * algorithmic parameter, not a memory knob) and the R* ratios that still feed
 * NUMSITES/NUMEXONS/MAXBACKUP* in SetRatios (used for the dumpster hash size,
 * the "backup active" flag, and beggar.c's -B estimate).
 * ---------------------------------------------------------------------- */

/* Length of every processed fragment       */
#define LENGTHSi 500000

/* Overlap between 2 fragments              */
#define OVERLAP 10000

/* One signal per L / RSITES bp             */
#define RSITES 1

/* /\* One U12 signal per L / RU12SITES bp             *\/ */
/* #define RU12SITES 6   */

/* One exon per L / REXONS bp  (divisor: smaller -> more exon headroom) */
#define REXONS 0.75

/* /\* One U12 intron-flanking exon per L / RU12EXONS bp               *\/ */
/* #define RU12EXONS 6  */

/* Estimated amount of backup signals (divisor: smaller -> more backup) */
#define RBSITES 75

/* Estimated amount of backup exons   (divisor: smaller -> more backup) */
#define RBEXONS 125

/* Ratios for every exon type               */
#define RFIRST 3
#define RINTER 1
#define RTERMI 2
#define RSINGL 4
#define RORF   4
#define RUTR   0.5

/* Exon-sort table factor (only used now by beggar.c's -B estimate) */
#define FSORT 12

/* Basic values (in addition to ratios)     */
#define BASEVALUESITES_SHORT 100000
#define BASEVALUEEXONS_SHORT 6000
#define BASEVALUESITES_LARGE 600000
#define BASEVALUEEXONS_LARGE 600000

/* Max number of annotations per locus      */
#define MAXEVIDENCES 200000       
#define MAXSITESEVIDENCES 3*MAXEVIDENCES

/* Max number of HSP per locus/frame/strand */
#define MAXHSP 40000000             

/* UTR prediction parameters */            
#define MAXUTRDONORS 8
#define MAX3UTREXONLENGTH 5000
#define MAXUTREXONLENGTH 1000
/* UTRMAXGAP used to be 55 */ 
#define UTRMAXGAP 5

/* read params */
#define RREADS 1
#define COV 15

/* Length of allowed UTR including stop codon before intron: used to be 55 */
#define MAXNMDLENGTH 1000

/* Max number of locus in multi-fasta files */
#define MAXNSEQUENCES 100         

/* Initial capacity of the growable per-gene info[] array (CookingGenes).   */
/* Not a ceiling: the array doubles on demand, so gene count is data-driven. */
#define INITGENES 256

/* Per-exon translation guard: max amino acids one exon contributes in       */
/* Translate() (a single exon is always far below this; PrintExons relies on  */
/* it too). The WHOLE-protein length is no longer capped -- see INITAA.       */
#define MAXAA 50000

/* Initial capacity of the growable whole-protein buffer (prot/sAux in        */
/* TranslateGene). Not a ceiling: it grows to fit the assembled protein.      */
#define INITAA 8192

/* Initial capacity of the growable cDNA/tDNA transcript buffers (GetcDNA/   */
/* GetTDNA). Not a ceiling: the buffers grow to fit the transcript length.   */
#define INITCDNA 4096

/* Initial capacity of the growable exon-sort table (SortExons). It grows on  */
/* demand, so the table is no longer reserved at FSORT*NUMEXONS up front nor  */
/* capped (the old "increase FSORT" hard error is gone).                      */
#define INITSORT 4096

/* Initial capacity of the growable per-type exon Build arrays (packExons).   */
/* Not a ceiling: each Build* grows its array on demand, so the old           */
/* "decrease RFIRST/RINTER/..." hard errors and the NUMEXONS/Rxxx reservation */
/* are gone.                                                                  */
#define INITEXONS 4096

/* Initial capacity of the growable site-sort scratch buffers (SortSites,     */
/* via packSortSites). Not a ceiling: each grows on demand, so the old        */
/* "increase FSORT" abort on the site-sort path is gone.                      */
#define INITSITESORT 1024

/* Initial capacity of the growable predicted-site arrays (packSites). Not a   */
/* ceiling: the site builders grow each array on demand, so the old NUMSITES  */
/* reservation and the "decrease RSITES" hard errors are gone.                */
#define INITSITES 4096

/* Initial capacity of each growable sort-by-donor array (packGenes->d[class], */
/* filled by BuildSort). Grows on demand, replacing the FDARRAY*NUMEXONS       */
/* reservation (and an unchecked overflow).                                   */
#define INITDARRAY 256

/* Chunk size (entries) of the growable dumpster backup arrays (BackupGenes.c).*/
/* Chunks are stable-address: allocated once, never moved, so backed-up        */
/* exons/sites keep their pointers valid while the dumpster grows on demand.   */
#define DUMPCHUNK 4096

/* Maximum number of isochores              */
#define MAXISOCHORES 4           

/* Minimum exon length (internal exons)     */
#define EXONLENGTH 18            

/* Minimum single gene length               */
#define SINGLEGENELENGTH 60      

/* Minimum ORF length                       */
#define ORFLENGTH 60             

/* Min length to compute coding potential   */
#define MINEXONLENGTH 6          

/* Score for exons with L <  MINEXONLENGTH  */
#define MINSCORELENGTH 0.0       

/* Region around exon to measure G+C        */
#define ISOCONTEXT 1000          

/* Number of nucloetides to scan for PPTs   */
/* or Branch Points before the Acceptor site*/
/* #define ACCEPTOR_CONTEXT 25*/
#define ACCEPTOR_CONTEXT 50
#define PPT_ACC_MAXDIST 14

/* Minimum distance between branch point    */
/* and acceptor site                        */
#define MIN_U12BPACC_DIST 7
#define MIN_U2BPACC_DIST 15
#define OPT_U12BP_DIST 12
#define OPT_U2BP_DIST 25

#define U12BP_PENALTY_SCALING_FACTOR 6 /* used to be 15 */
#define U2BP_PENALTY_SCALING_FACTOR 0 /* used to be 15 */

/* Recursive splice site thresholds */
#define RDT 4 /*donor*/
#define RAT 4 /*acceptor*/
#define sRSSMARKOVSCORE "RSS_Markov_Score"
#define sRSS_DONOR_SCORE_CUTOFF "RSS_Donor_Score_Cutoff"
#define sRSS_ACCEPTOR_SCORE_CUTOFF "RSS_Acceptor_Score_Cutoff"

/* Markov score penalty for unknown symbols */
#define NULL_OLIGO_SCORE  -4     

/* Maximum profile dimension (length)       */
#define PROFILEDIM 100           

/* Maximum number of chars (locus names)    */
#define LOCUSLENGTH 500          

/* Maximum oligo (word) length (Markov)     */
#define OLIGOLENGTH 10           

/* Maximum chars per input line             */
#define MAXLINE 1000             

/* Characters per fasta line                */
#define FASTALINE 60             

/* Maximum length for strings (mess)        */
#define MAXSTRING 600            

/* Mark rules up as blocking in Gene model  */
#define BLOCK 1                  
#define NONBLOCK 0               

/* Array range in C: 0..N-1                 */
#define COFFSET 1                

/* Dictionary definitions (hash)            */
#define MAXENTRY 97              
#define MAXTYPE 50               
#define MAXINFO 100
#define NOTFOUND -1
#define UNKNOWNAA 'X'

/* Dumpster hash table size (factor)        */ 
#define HASHFACTOR 3             

/* maximum length of filenames              */
#define FILENAMELENGTH 5000       

/* Name of default parameter file           */
#define PARAMETERFILE  "param.default"   

/* Constants:                               */
#define FRAMES 3                   
#define STRANDS 2                      
#define LENGTHCODON 3                  
#define PERCENT 100
#define MINUTE 60                    
#define MEGABYTE 1048576
#define MAXTIMES 100
#define PROT 0
#define TDNA 2
#define DNA  1

/* Strands                                  */
#define cFORWARD '+'
#define cREVERSE '-'

#define FORWARD 0                      
#define REVERSE 1

#define sFORWARD "Forward"
#define sREVERSE "Reverse"

#define xmlFORWARD "fwd"
#define xmlREVERSE "rvs"

/* Signals                                  */
#define ACC 0                    
#define DON 1
#define STA 2
#define STO 3
#define POL 4
#define TSS 5
#define TES 6

#define sACC "Acceptor"
#define sDON "Donor"
#define sSTA "Start"
#define sSTO "Stop"
#define sPOL "PolyA"
#define sTSS "TSS"
#define sTES "TES"
#define sPPT "PolyPyrimidineTract"
#define sBP  "BranchPoint"

/* Splice Classes                                  */
#define U2 0                    
#define U12gtag 1
#define U12atac 2

/* Intron Types                           */
#define sU2type "U2"
#define sU12type "U12"

/* Intron Subtypes                           */
#define MAXSUBTYPE 10
#define MAXSPLICETYPE 5
#define sU2 "U2"  
#define sU2gcag "U2gcag"            
#define sU12gtag "U12gtag"
#define sU12atac "U12atac"
#define sU2gta "U2gta"
#define sU2gtg "U2gtg"
#define sU2gty "U2gty"

/* Header profiles                          */
#define sMarkov "Markov_order"
#define sprofilePolyA "PolyA_Signal_profile"
#define sprofilePPT "Poly_Pyrimidine_Tract_profile"
#define sprofileBP  "Branch_point_profile"
#define sprofileACC "Acceptor_profile"
#define sprofileU12BP  "U12_Branch_point_profile"
#define sprofileU12gtagACC "U12gtag_Acceptor_profile"
#define sprofileU12atacACC "U12atac_Acceptor_profile"
#define sprofileDON "Donor_profile"
#define sprofileU2gcagDON "U2gcag_Donor_profile"
#define sprofileU12gtagDON "U12gtag_Donor_profile"
#define sprofileU12atacDON "U12atac_Donor_profile"
#define sprofileU2gtaDON "U2gta_Donor_profile"
#define sprofileU2gtgDON "U2gtg_Donor_profile"
#define sprofileU2gtyDON "U2gty_Donor_profile"

/* U12 splice site (sum d + a) and exon (sum d + a) score thresholds    */
#define sU12_SPLICE_SCORE_THRESH "U12_Splice_Score_Threshold"
#define sU12_EXON_SCORE_THRESH "U12_Exon_Score_Threshold"
/* U12 acceptor splice site initial threshold factor (added to final threshold given in param file) */
#define U12ACC_CUTOFF_FACTOR -7

/* Header exon weights                          */
#define sExon_weights "Exon_weights"
#define sU12_EXON_WEIGHT "U12_Exon_weight"

/* Header evidence factor and weight */
#define sEVIDENCEF "Evidence_Factor"
#define sEVIDENCEW "Evidence_Exon_Weight"
#define sBKGD_SUBTRACT_FLANK_LENGTH "BKGD_SUBTRACT_FLANK_LENGTH"
#define sNUMISO "number_of_isochores"

/* Exons                                    */
#define FIRST             0
#define INTERNAL          1
#define TERMINAL          2
#define SINGLE            3
#define ORF               4
#define ZEROLENGTH        5
#define INTRON            6
#define UTRFIRST          7
#define UTRFIRSTHALF      8
#define UTRINTERNAL       9
#define UTR5INTERNALHALF 10
#define UTR3INTERNALHALF 11
#define UTRTERMINALHALF  12
#define UTRTERMINAL      13
#define UTRINTRON        14
#define UTR5INTRON       15
#define UTR3INTRON       16

#define sFIRST    "First"              
#define sINTERNAL "Internal"
#define sINTRON   "Intron"
#define sZEROLENGTH "RSS"
#define sTERMINAL "Terminal"
#define sSINGLE   "Single"
#define sORF      "ORF"
#define sEXON     "Exon"
#define sPROMOTER "Promoter"              
#define sUTRFIRST          "UTR_First"
#define sUTRFIRSTHALF      "UTR_First_Half"
#define sUTRINTERNAL       "UTR_Internal"
#define sUTR5INTERNALHALF  "UTR_5prime_Internal_Half"
#define sUTR3INTERNALHALF  "UTR_3prime_Internal_Half"
#define sUTRTERMINALHALF   "UTR_Terminal_Half"
#define sUTRTERMINAL       "UTR_Terminal"
#define sUTRINTRON         "UTR_Intron"
#define sUTR5INTRON        "UTR_5prime_Intron"
#define sUTR3INTRON        "UTR_3prime_Intron"

#define sBEGIN    "Begin"
#define sBEGINFWD "Begin+"
#define sBEGINRVS "Begin-"
#define sEND      "End"
#define sENDFWD   "End+"
#define sENDRVS   "End-"

/* Infinity: positions in sequence          */
#define INFI 2147483647        

/* Infinity: word in Gene model             */
#define SINFI "Infinity"         

/* Infinity: score functions                */
#define INF  1.7976931348623157E+308  

/* Score (annotations with score=".")       */
#define MAXSCORE 10000.0         

/* ReadSeq: Message (info)/X chars          */
#define MESSAGE_FREQ 100000      

/* Field group in gff: Not grouped exons    */
#define NOGROUP "NON_GROUPED"    

#define NOVALUE 0              

/* Estimation values: memory account        */
#define AVG_DIM   15
#define AVG_ORDER  3
#define OLIGO_DIM  5

/* Macros (functions)                       */
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))


/*************************************************************************
B. DATA TYPES
*************************************************************************/

typedef struct s_node *pnode;          
typedef struct s_node
{
  char s[MAXSTRING];
  int key;
  pnode next;
} node; 

typedef struct s_dict                  
{
  pnode T[MAXENTRY];
  int nextFree; 
} dict;

typedef struct s_profile               
{
  int    dimension;
  int    offset;
  float  cutoff;
  int    order;
  float  afactor;
  float  bfactor;
  int    acc_context;
  int    dist;
  int    opt_dist;
  float    penalty_factor;

  long dimensionTrans;
  float*  transitionValues[PROFILEDIM];
} profile;

/* A single predicted splice site / start / stop / TSS / TES signal. */
typedef struct s_site
{
  long Position;         /* coordinate in the (fragment-local) sequence */
  float Score;            /* final signal score (after all sub-profiles below) */
  float ScoreAccProfile;  /* acceptor PWM/Markov score alone, before adding BP/PPT */
  float ScoreBP;          /* branch-point sub-profile score (acceptors only) */
  float ScorePPT;         /* poly-pyrimidine-tract sub-profile score (acceptors only) */
  int PositionBP;         /* branch point position, relative to this Acceptor */
  int PositionPPT;        /* poly-pyrimidine-tract position, relative to this Acceptor */
  char type[MAXSPLICETYPE];    /* e.g. "Acceptor", "Donor", "Start", "Stop" (site kind) */
  char subtype[MAXSUBTYPE];    /* e.g. "U2", "U12gtag", "U12atac", "GC-AG" (signal flavor) */
  short class;            /* U2 / U12gtag / U12atac (see #define's above); drives Ga's 3rd axis */
} site;

/* All predicted signals for one strand of one fragment, one array per
   signal kind (TS/TE = transcription start/end sites, UTR-only). */
typedef struct s_packSites
{
  site* StartCodons;
  site* AcceptorSites;
  site* DonorSites;
  site* StopCodons;
  site* TS;
  site* TE;

  long  nStartCodons;
  long  nAcceptorSites;
  long  nDonorSites;
  long  nStopCodons;
  long  nTS;
  long  nTE;

  /* Allocated capacity of each growable site array (grown by the Build/Get
     site builders; the arrays are no longer reserved at NUMSITES nor capped). */
  long  capStartCodons;
  long  capAcceptorSites;
  long  capDonorSites;
  long  capStopCodons;
  long  capTS;
  long  capTE;

  long nSites;
} packSites;

/* Scratch buffers SortSites uses to stage sorted sites before copying them
   back into the originating packSites arrays (donor/acceptor/TS/TE). Each
   grows on demand (see GrowSiteArray); capacity persists across fragments. */
typedef struct s_packSortSites
{
  site* donorsites;
  long  donorsitescap;
  site* acceptorsites;
  long  acceptorsitescap;
  site* tssites;
  long  tssitescap;
  site* tesites;
  long  tesitescap;
} packSortSites;

/* A single predicted (or annotated/evidence) exon, and also the node type
   used to chain exons into genes: genamic's DP walks PreviousExon links to
   build up a gene, scoring each addition into GeneScore (see genamic.c). */
typedef struct s_exonGFF *pexonGFF;
typedef struct s_exonGFF
{
  site* Acceptor;         /* left (5') splice/start signal bounding this exon */
  site* Donor;            /* right (3') splice/stop signal bounding this exon */
  char Type[MAXTYPE];     /* "First"/"Internal"/"Terminal"/"Single"/"Intron"/UTR variants/... */
  short Frame;            /* reading frame this exon starts in (0/1/2) */
  short Remainder;        /* bases left over from an incomplete trailing codon */
  char Strand;            /* '+' / '-' / '*' (the Ghost/placeholder exon) */
  float PartialScore;     /* the coding-potential (Markov) component of Score */
  float HSPScore;         /* the homology/HSP-support component of Score */
  float R;                /* RNA-seq coverage evidence (e.g. RPKM) for this exon, if any */
  float Score;             /* this exon's own score (site + coding + HSP + weight, see ScoreExons) */
  pexonGFF PreviousExon;  /* best predecessor exon in the assembled gene chain (Ghost if none) */
  double GeneScore;       /* PreviousExon->GeneScore + Score: best gene score ending at this exon */
  char Group[MAXSTRING];  /* evidence/annotation group name (NOGROUP if none); see genamic's block rule */
  int offset1;            /* Acceptor-position correction (fragment-splitting/evidence coordinate offset) */
  int offset2;            /* Donor-position correction (fragment-splitting/evidence coordinate offset) */
  short lValue;           /* left split-codon signature (set by ComputeStopInfo; see genamic.c) */
  short rValue;           /* right split-codon signature (set by ComputeStopInfo; see genamic.c) */
  short evidence;         /* 1 if this exon comes from external evidence/annotation, 0 if predicted */
  short selected;         /* scratch mark used by SwitchFrames to track already-processed d-array entries */
  short three_prime_partial; /* 1 if the gene's 3' end is cut off by the sequence/fragment boundary (GFF3 output) */
  short five_prime_partial;  /* 1 if the gene's 5' end is cut off by the sequence/fragment boundary (GFF3 output) */
} exonGFF;

/* One array per exon type, all built and filled per fragment/strand by the
   Build*Exons functions (see e.g. BuildInitialExons.c); UTR arrays only
   used when UTR prediction (-u) is on. Each array grows on demand (see the
   matching cap* field and GrowExonArray). */
typedef struct s_packExons
{
  exonGFF* InitialExons;
  exonGFF* InternalExons;
  exonGFF* ZeroLengthExons;
  exonGFF* TerminalExons;
  exonGFF* Singles;
  exonGFF* ORFs;
  exonGFF* UtrInitialExons;
  exonGFF* UtrInitialHalfExons;
  exonGFF* UtrInternalExons;
  exonGFF* Utr5InternalHalfExons;
  exonGFF* Utr3InternalHalfExons;
  exonGFF* UtrTerminalHalfExons;
  exonGFF* UtrTerminalExons;
  long nInitialExons;                  
  long nInternalExons;
  long nZeroLengthExons;
  long nTerminalExons;
  long nSingles;
  long nORFs;
  long nUtrInitialExons;
  long nUtrInitialHalfExons;
  long nUtrInternalExons;
  long nUtr5InternalHalfExons;
  long nUtr3InternalHalfExons;
  long nUtrTerminalHalfExons;
  long nUtrTerminalExons;
  long nExons;
  /* Allocated capacity of each growable per-type array (grown by Build*) */
  long capInitialExons;
  long capInternalExons;
  long capZeroLengthExons;
  long capTerminalExons;
  long capSingles;
  long capORFs;
  long capUtrInitialExons;
  long capUtrInitialHalfExons;
  long capUtrInternalExons;
  long capUtr5InternalHalfExons;
  long capUtr3InternalHalfExons;
  long capUtrTerminalHalfExons;
  long capUtrTerminalExons;
} packExons;

/* The gene-assembly DP state, persisted across all fragments of one locus
   (reset only when moving to the next input sequence). See the header
   comment of genamic.c for how these are used together. */
typedef struct s_packGenes
{
  exonGFF* ***Ga;   /* Ga[class][frame][spliceclass]: best partial gene ending in that DP cell */
  exonGFF* Ghost;   /* placeholder "empty gene" every Ga cell starts pointing at (Strand '*') */
  exonGFF* GOptim;  /* best-scoring complete gene found so far, across the whole locus */
  exonGFF* **d;     /* d[class]: exons compatible with gene-model class `class`, sorted by Donor position */
  long* dcap;       /* allocated capacity of each d[class] (grown by BuildSort) */
  long* km;         /* km[class]: number of exons currently in d[class] */
  long* je;         /* je[class]: how far d[class] has been folded into Ga so far (forward-only cursor) */
} packGenes;

typedef struct s_packEvidence
{
  site* vSites;
  long nvSites;
  exonGFF* vExons; 
  long nvExons;  
} packEvidence;

typedef struct s_HSP
{
  long Pos1;
  long Pos2;
  float Score;
} HSP;

typedef struct s_packHSP
{
  HSP*** sPairs;
  long* nSegments;
  long nTotalSegments;
  int visited;
} packHSP;

typedef struct s_packExternalInformation
{
  dict* locusNames;

  packHSP** homology;
  packEvidence** evidence;

  long nSequences;
  long nvExons;
  long nHSPs;

  long i1vExons;
  long i2vExons;
  long ivExons;
  
  long* iSegments;
  float** sr;
  float** readcount;
} packExternalInformation;

/* Hash-bucket entry identifying one already-backed-up exon (see DumpHash.c);
   the key fields mirror exonGFF's Acceptor/Donor/Frame/Strand/Type so
   getExonDumpHash can recognize a duplicate without re-hashing the exon
   from scratch. asub/dsub are unused (their fill-in is commented out in
   setExonDumpHash) -- kept only for the on-disk struct layout. */
typedef struct s_dumpNode *pdumpNode;
typedef struct s_dumpNode
{
  long acceptor;    /* the exon's Acceptor->Position */
  long donor;       /* the exon's Donor->Position */
  short aclass;     /* the exon's Acceptor->class (U2/U12gtag/U12atac) */
  short dclass;     /* the exon's Donor->class */
  char asub[MAXSUBTYPE];
  char dsub[MAXSUBTYPE];
  short frame;      /* the exon's Frame */
  char strand;      /* the exon's Strand */
  char type[MAXTYPE]; /* the exon's Type */

  exonGFF* exon;    /* the actual backed-up exon (a stable pointer into a dumpster chunk) */
  pdumpNode next;   /* next node in this hash bucket's collision chain */
} dumpNode;

/* Hash table (fixed MAXDUMPHASH buckets, chained) so backupGene can tell in
   O(1) whether an exon has already been backed up into the dumpster, instead
   of rescanning it. */
typedef struct s_dumpHash
{
  pdumpNode* T;   /* T[0..MAXDUMPHASH): bucket array, each a collision chain */
  long total;     /* total number of exons currently hashed (diagnostic only) */
} dumpHash;

typedef struct s_packDump
{
  /* Chunked, stable-address backup arrays (see BackupGenes.c). dumpSites and
     dumpExons are arrays of fixed-size chunks; a chunk is never moved once
     allocated, so a backed-up site/exon keeps its address for the whole run
     while the arrays grow on demand -- replacing the old fixed ring buffer
     that recycled (and could overwrite still-referenced backups). */
  site** dumpSites;
  long ndumpSites;
  long dumpSitesChunks;

  exonGFF** dumpExons;
  long ndumpExons;
  long dumpExonsChunks;

  dumpHash* h;
} packDump;

typedef struct s_account               
{
  long starts, starts_r,
    stops, stops_r,
    acc, acc_r,
    don, don_r,
    tss, tss_r,
    tes, tes_r;

  long first, first_r,
    internal, internal_r,
    terminal, terminal_r,
    single, single_r,
    orf, orf_r, zle, zle_r, utr, utr_r;

  long totalExons;
  

  int tSites, tExons, tGenes, tSort, tScore, tBackup;
  time_t tStart;

} account;

typedef struct s_packGC
{
  long* GC;
  long* N;
} packGC;

typedef struct s_paramexons            
{
  float siteFactor;

  float exonFactor;
  float OligoCutoff;

  float HSPFactor;

  float ExonWeight;
/*   float U12atacExonWeight; */
/*   float U12gtagExonWeight; */
  float ExonCutoff;
} paramexons;

typedef struct s_gparam                
{
  int leftValue;
  int rightValue;

  profile* PolyASignalProfile;               
  profile* StartProfile;               
  profile* AcceptorProfile;
  profile* PolyPTractProfile;
  profile* BranchPointProfile;
  profile* DonorProfile;
  profile* U2gcagDonorProfile;
  profile* U2gtaDonorProfile;
  profile* U2gtgDonorProfile;
  profile* U2gtyDonorProfile;
  profile* U12gtagAcceptorProfile;
  profile* U12BranchPointProfile;
  profile* U12gtagDonorProfile;
  profile* U12atacAcceptorProfile;
  profile* U12atacDonorProfile;
  profile* StopProfile;

  float* OligoLogsIni[3];              
  float* OligoLogsTran[3];

  long OligoDim;                       
  long OligoDim_1;           
  int  OligoLength;

  paramexons* Initial;                 
  paramexons* Internal;
  paramexons* Terminal;
  paramexons* Single;
  paramexons* utr;
                                                                            
  float* OligoDistIni[FRAMES]; 
  float* OligoDistTran[FRAMES];
  
  int MaxDonors;                       

  dict* D;                             
  int*  nc;
  int*  ne;
  long* md;
  long* Md;
  int UC[MAXENTRY][MAXENTRY];
  int DE[MAXENTRY][MAXENTRY];
  int block[MAXENTRY];
  int nclass;

} gparam;


/*************************************************************************
C. IMPORTED HEADERS
*************************************************************************/
float strtof(const char *nptr, char **endptr);

void PrintExonGFF(exonGFF *e, char* Name, char* Source);

void PrintGeneGFF(exonGFF *e, char* Name, char* Source);

void printError(char *s);

void printMess(char* s);

void printRes(char* s);

void printReadingInfo(char* s);

long GetSitesWithProfile(char *s, profile *p, site **st, long* cap, long l1, long l2);

long GetTSS(
	    site** sc, long* cap,
	    site* Acceptors, long nAcceptors,
	    packExternalInformation* external,
	    packHSP* hsp,
	    int Strand,
	    long LengthSequence,
	    long l1,
	    long l2
	    );

long GetTES(
	    site** sc, long* cap,
	    site* Donors, long nDonors,
	    packExternalInformation* external,
	    packHSP* hsp,
	    int Strand,
	    long LengthSequence,
	    long l1,
	    long l2,
	    long ns
	    );

long BuildDonors(char* s,short class,char* type,
		 char* subtype,
		 profile* p,
		 site** st,
		 long l1,
		 long l2,
		 long ns,
		 long* cap,
		 int Strand,
		 packExternalInformation* external
		 );
float PeakEdgeScore(long Position, 
		    int Strand, 
		    packExternalInformation* external, 
		    long l1, long l2, int win);
int ClusterEdge(long Position, 
		int Strand, 
		packExternalInformation* external, 
		long l1, long l2);
long GetStopCodons(char *s, profile *p, site **sc, long* cap, long l1, long l2);

long BuildInitialExons(site *Start, long nStarts, 
                       site *Donor, long nDonors,
                       site *Stop, long nStops,
                       int MaxDonors,
		       char* ExonType,
		       char* Sequence,
                       exonGFF **Exon, long* cap );

long BuildInternalExons(site *Acceptor, long nAcceptors,
                        site *Donor, long nDonors,
                        site *Stop, long nStops,
                        int MaxDonors,
			char* ExonType,
			char* Sequence,
                        exonGFF** Exon, long* cap);

long BuildZeroLengthExons(site *Acceptor, long nAcceptors,
                        site *Donor, long nDonors,
                        site *Stop, long nStops,
                        int MaxDonors,
			char* ExonType,
			char* Sequence,
                        exonGFF** Exon, long* cap);

long BuildTerminalExons (site *Acceptor, long nAcceptors,
                         site *Stop, long nStops,
                         long LengthSequence,
                         long cutPoint,
			 char* ExonType,
			 char* Sequence,
			 exonGFF** Exon, long* cap);

long BuildSingles(site *Start, long nStarts,
                  site *Stop, long nStops,
                  long cutPoint,
		  char* Sequence,
                  exonGFF **Exon, long* cap);

long BuildORFs(site *Start, long nStarts,
	       site *Stop, long nStops,
	       long cutPoint,
	       char* Sequence,
	       exonGFF **Exon, long* cap);

long BuildUTRExons(
		   site *Start, long nStarts,
		   site *Donor, long nDonors,
		   int MaxDonors,
		   int MaxExonLength,long cutPoint,
		   char* ExonType,
		   exonGFF **Exon, long* cap);

packSites* RequestMemorySites();
packExons* RequestMemoryExons();
exonGFF* RequestMemorySortExons();

/* Grow an exonGFF array *buf to hold at least `need` entries (realloc-double,
   zeroing the new entries since the arrays were originally calloc'd). Shared
   by the growable per-type Build arrays and the exon-sort table. */
void GrowExonArray(exonGFF** buf, long* cap, long need);
site* RequestMemorySortSites();

/* Grow a site array *buf to hold at least `need` entries (realloc-double,
   zeroing the new entries). Used by the growable site-sort scratch buffers. */
void GrowSiteArray(site** buf, long* cap, long need);
gparam* RequestMemoryParams();
packGenes* RequestMemoryGenes();
packDump* RequestMemoryDumpster();
dict* RequestMemoryAaDictionary();
void RequestMemoryProfile(profile* p);
account* RequestMemoryAccounting();
packExternalInformation* RequestMemoryExternalInformation();

void readargv (int argc,char *argv[],
	       char *ParamFile, char* SequenceFile,
	       char *ExonsFile, char* HSPFile, char* GenePrefix);

int readparam (char *name, gparam** isochores);

account* InitAcc();

void OutputHeader(char* locus, long l);

int IniReadSequence(FILE* seqfile, char* line);

int ReadSequence (FILE* seqfile, char* Sequence, char* nextLocus);

long FetchSequence(char *s, char* r);

long ReadExonsGFF (char *FileName, 
		   packExternalInformation* external, 
		   dict* d);

void SwitchPositions(packExons *allExons);

void SearchEvidenceExons(packExternalInformation* external,
			 packEvidence* pv, 
			 long l2);

void SortExons(packExons* allExons,
               packExons* allExons_r,
	       packExternalInformation* external,
	       packEvidence* pv,
               exonGFF** Exons,
               long* Exonscap,
               long l1, long l2,long lowerlimit,
	       long upperlimit);

void SortSites(site *Sites, long nSites, site **sortedSites, long* cap,
               long l1, long l2
	       );

void SwitchCounters(packExternalInformation* external);

void Output(packSites* allSites, packSites* allSites_r,
            packExons* allExons, packExons* allExons_r,
            exonGFF* exons, long nExons, char* Locus, 
	    long l1, long l2, long lowerlimit, char* Sequence, gparam* gp, dict* dAA, char* GenePrefix);

void updateTotals(account *m,
                  packSites* allSites,
                  packSites* allSites_r,
                  packExons* allExons,
                  packExons* allExons_r);

void genamic(exonGFF *E, long nExons, packGenes* pg, gparam* gp);

void BackupGenes(packGenes* pg, int nclass, packDump* d);

void BackupArrayD(packGenes* pg, long accSearch,
                  gparam* gp, packDump* dumpster);

void cleanGenes(packGenes* pg, int nclass, packDump* dumpster);

void cleanDumpHash(dumpHash *h);

void OutputGene(packGenes* pg, long nExons, char* Locus,
                char* Sequence, gparam* gp, dict* dAA, char* GenePrefix);

void OutputStats(char* Locus);
void OutputTime();

void RecomputePositions(packSites* allSites, long l);

void cleanAcc(account* m);

void PrintProfile (profile *p, char* signal);

long ReadGeneModel (FILE *file, 
		    dict *d, 
		    int nc[], 
		    int ne[], 
                    int UC[][MAXENTRY], 
		    int DE[][MAXENTRY], 
		    long md[], 
		    long Md[], 
		    int block[]);

long ForceGeneModel (dict* d,
		     int nc[], int ne[],
		     int UC[][MAXENTRY],
		     int DE[][MAXENTRY],
		     long md[], long Md[],
		     int block[]);

void PrintSites (site *s, long ns,int type,
                 char Name[], int Strand,
                 long l1, long l2, long lowerlimit,
                 char* seq,
                 profile *p);

void PrintExons (exonGFF *e, long ne, int type, char Name[],
                 long l1, long l2, char* Sequence,  dict* dAA, char* GenePrefix); 

void resetDict(dict* d);

int setkeyDict(dict *d, char s[]);

int getkeyDict(dict *d, char s[]);

int Translate(long p1, long p2, short fra, short rmd,
              char *s, dict* dAA, char sAux[]);

void ReverseSubSequence(long p1, long p2, char* s, char* r);

void CorrectExon(exonGFF *e);

void CorrectUTR(exonGFF *e);

void CorrectORF(exonGFF *e);

void SwitchFrames(exonGFF* e, long n);

void SwitchFramesDa(packGenes* pg, int nclass);

void SwitchFramesDb(packGenes* pg, int nclass);

void UndoFrames(exonGFF *e, long n);

void BuildSort(dict *D, int nc[], int ne[], int UC[][MAXENTRY],
               int DE[][MAXENTRY], int nclass, long km[], long dcap[],
	       exonGFF* **d, exonGFF *E, long nexons);

void PrintSite(site *s, int type, char Name[], int Strand,
               char* seq, profile *p);

void PrintGCDS(exonGFF *e, char Name[], char* s, dict* dAA, 
		long ngen, int AA1, int AA2, int nAA,
		int numInt, char* GenePrefix);

void PrintGUTR(exonGFF *e, char Name[], char* s, long ngen,
		int numInt, char* GenePrefix);

void PrintGExon(exonGFF *a, int nSegments, char Name[],long ngen,
		  int numInt, char* GenePrefix, int evidence, float score);
			   
void PrintGIntron(exonGFF *d, exonGFF *a, char Name[],long ngen,
		  int numInt, char* GenePrefix, int evidence, float score, char* eType);

void PrintGGene(exonGFF *s, exonGFF *e, char Name[],
		long ngen, float score, char* GenePrefix);

void PrintGmRNA(exonGFF *s, exonGFF *e, char Name[],
		long ngen, float score, char* GenePrefix);


void PrintXMLExon(exonGFF *e, char Name[], 
		  long ngen, long nExon, 
		  int type1, int type2, 
		  char* GenePrefix);

void TranslateGene(exonGFF* e,
                   char* s,
                   dict* dAA,
                   long nExons,
                   int** tAA,
                   char** prot,
                   long* cap,
                   long* nAA);

void resetDumpHash(dumpHash* h);

void setExonDumpHash(exonGFF* E, dumpHash *h);

exonGFF* getExonDumpHash(exonGFF* E, dumpHash *h);
 
void setAADict(dict *d, char s[], char aA);

char getAADict(dict *d, char s[]);

void CookingGenes(exonGFF *e, char Name[], char* s,
                  gparam* gp, dict* dAA, char* GenePrefix);

float MeasureSequence(long l1,long l2,char* s);

gparam ** RequestMemoryIsochoresParams();

long ReadHSP (char* FileName, packExternalInformation* external);

void shareGeneModel(gparam** isochores, int n);

long OligoToInt(char *s, int ls, int cardinal);

char* RequestMemorySequence(long L);

void ProcessHSPs(long l1,
                long l2,
                int Strand,
		packExternalInformation* external,
                packHSP* hsp);

void ScoreExons(char *Sequence, 
                packExons* allExons, 
                long l1,
                long l2,
                int Strand,
		packExternalInformation* external,
                packHSP* hsp,
                gparam** isochores,
                int nIsochores,
                packGC* GCInfo);

void GetcDNA(exonGFF* e, char* s, long nExons, char** cDNA, long* cap, long* nNN);
void GetTDNA(exonGFF* e, char* s, long nExons, char** tDNA, long* cap, long* nTN);

float ComputeGC(packGC* GCInfo, long inigc, long endgc);

void GCScan(char* s, packGC* GCInfo, long l1, long l2);

void beggar(long L);

void SetRatios(long* NUMSITES,
               long* NUMEXONS,
               long* MAXBACKUPSITES,
               long* MAXBACKUPEXONS,
               long L);

long analizeFile(char* SequenceFile);

packGC* RequestMemoryGC();

int SelectIsochore(float percent, gparam** isochores);

void  manager(char *Sequence, long LengthSequence,
	      packSites* allSites,
	      packExons* allExons,
	      long l1, long l2,long lower, long upper,
	      int Strand,
	      packExternalInformation* external,
	      packHSP* hsp,
	      gparam* gp,
	      gparam** isochores,
	      int nIsochores,
	      packGC* GCInfo,
	      packSortSites* sortSites
	      );

void resetEvidenceCounters(packExternalInformation* external);

void ComputeStopInfo(exonGFF* e, char* s);

packHSP* SelectHSP(packExternalInformation* external,
		   char* Locus, 
		   long LengthSequence);

packEvidence* SelectEvidence(packExternalInformation* external,
			     char* Locus);

void SortHSPs(packHSP* p);

HSP* RequestNewHSP();

long  BuildU12Acceptors(char* s,
			short class,
			char* type,
			char* subtype,
			profile* u12_p,
			profile* u12bp,
			profile* ppt,
			site** st,
			long l1,
			long l2,
			long ns,
			long* cap,
			int Strand,
			packExternalInformation* external);

long  BuildAcceptors(char* s,
		     short class,
		     char* type,
		     char* subtype,
		     profile* p,
		     profile* ppt,
		     profile* bp,
		     site** st,
		     long l1,
		     long l2,
		     long ns,
		     long* cap,
		     int Strand,
		     packExternalInformation* external);
