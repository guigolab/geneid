/*************************************************************************
*                                                                        *
*   Module: geneid.h                                                     *
*                                                                        *
*   Definitions, data structures types and imported headers              *
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

/* $Id: geneid.h,v 1.12 2002-03-20 10:18:34 eblanco Exp $ */

/* Required libraries */
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <malloc.h>  
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <sys/stat.h>

/*************************************************************************
(A). DEFINITIONS
*************************************************************************/

#define LENGTHSi 100000          /* Length of every processed fragment   */ 
#define OVERLAP 10000            /* Overlap between 2 fragments          */

#define RSITES 5                 /* One signal per L / RSITES bp         */
#define REXONS 3                 /* One exon per L / REXONS bp           */

#define RBSITES 75               /* Estimated amount of backup signals   */
#define RBEXONS 125              /* Estimated amount of backup exons     */

#define RFIRST 1                 /* Ratios according to every exon type  */
#define RINTER 0.25
#define RTERMI 1
#define RSINGL 3 
#define RORF   3

#define FSORT 8                  /* Total number ox exons/split (factor) */
#define FDARRAY 5                /* Number of exons to save every  split */ 

#define BASEVALUESITES_SHORT 10000
#define BASEVALUEEXONS_SHORT 5000
#define BASEVALUESITES_LARGE 600000
#define BASEVALUEEXONS_LARGE 300000

#define MAXEVIDENCES 10000       /* Maximum number of annotations        */
#define MAXSR 10000              /* Maximum number of SR (strand/frame)  */

#define MAXGENE 20000            /* Maximum number of predicted genes    */
#define MAXEXONGENE 1000         /* Maximum number of exons in a gene    */
#define MAXAA 300000             /* Maximum lenght of a protein          */
#define MAXCDNA MAXAA*3          /* Maximum length of (cDNA) in genes    */
#define MAXISOCHORES 4           /* Maximum number of isochores          */

#define EXONLENGTH 18            /* minimum exon length (internal)       */
#define SINGLEGENELENGTH 60      /* minimum single gene length           */
#define ORFLENGTH 60             /* minimum ORF length                   */

#define MINEXONLENGTH 6          /* Min. length to score an exon (Markov)*/
#define MINSCORELENGTH 0.0       /* Score for <  MINEXONLENGTH           */

#define ISOCONTEXT 1000          /* Region around exon to measure G+C    */

#define NULL_OLIGO_SCORE  -4     /* Markov score penalty: N's            */
#define PROFILEDIM 100           /* maximum profile's dimension          */

#define LOCUSLENGTH 100          /* maximum number of chars (locus)      */
#define OLIGOLENGTH 10           /* maximum oligo (word) length          */
#define MAXLINE 1000             /* Maximum chars per input line         */
#define FASTALINE 60             /* Characters per fasta line            */
#define MAXSTRING 100            /* Maximum length for strings (mess)    */

#define BLOCK 1                  /* Mark a rule as blocking              */
#define NONBLOCK 0               /* Mark a rule as non blocking          */

#define COFFSET 1                /* Array range in C: 0..N-1             */

#define MAXENTRY 97              /* Dictionary definitions (hash)        */
#define MAXTYPE 50               /* Maximum number of chars/exonTypes    */
#define MAXINFO 100
#define NOTFOUND -1

#define HASHFACTOR 3             /* Dumpster hash table size (factor)    */ 

#define FILENAMELENGTH 200       /* maximum length of filenames          */
#define PARAMETERFILE  "param.default"   /* default parameters file      */

#define VERSION   "geneid_v1.1"  /* The name of the game                 */
#define SITES     "geneid_v1.1"  /* gff fields: source and feature       */
#define EXONS     "geneid_v1.1"        
#define EVIDENCE  "evidence"           

#define FRAMES 3                 /* Constants:                           */
#define STRANDS 2                      
#define LENGTHCODON 3                  
#define PERCENT 100
#define MINUTE 60                    
#define MEGABYTE 1048576
#define MAXTIMES 100
#define PROT 0
#define DNA  1

#define cFORWARD '+'
#define cREVERSE '-'

#define FORWARD 0                      
#define REVERSE 1

#define sFORWARD "Forward"
#define sREVERSE "Reverse"

#define xmlFORWARD "fwd"
#define xmlREVERSE "rvs"

#define ACC 0                    /* Signals                              */
#define DON 1
#define STA 2
#define STO 3

#define sACC "Acceptor"                
#define sDON "Donor"
#define sSTA "Start"
#define sSTO "Stop"

#define FIRST    0               /* Exons                                */
#define INTERNAL 1
#define TERMINAL 2
#define SINGLE   3
#define ORF      4

#define sFIRST    "First"              
#define sINTERNAL "Internal"
#define sTERMINAL "Terminal"
#define sSINGLE   "Single"
#define sORF      "ORF"
#define sEXON     "Exon"

#define INFI 999999999           /* the biggest in geneid world          */
#define SINFI "Infinity"         
#define INF  1.7976931348623157E+308  

#define MAXSCORE 10000.0         /* Score (forced annotations)           */

#define MESSAGE_FREQ 100000      /* Message (info) per input X chars     */

#define NOGROUP "NON_GROUPED"    /* Field group in gff: Not grouped exons*/

#define NOVALUE 0              

/* Macros (functions) */
#define MIN(a,b) (a<b)?a:b;
#define MAX(a,b) (a>b)?a:b;


/*************************************************************************
(B). DATA TYPES
*************************************************************************/

typedef struct s_profile               
{
  int    dimension;
  int    offset;
  float  cutoff;
  int    order;

  long dimensionTrans;
  float*  transitionValues[PROFILEDIM];
} profile;

typedef struct s_site                  
{
  long Position;                       
  double Score;                        
} site;

typedef struct s_packSites             
{
  site* StartCodons;                   
  site* AcceptorSites;
  site* DonorSites;
  site* StopCodons;

  long  nStartCodons;                  
  long  nAcceptorSites;
  long  nDonorSites;
  long  nStopCodons;

  long nSites;                         
} packSites;

typedef struct s_exonGFF *pexonGFF;    
typedef struct s_exonGFF
{
  site* Acceptor;
  site* Donor;
  char Type[MAXTYPE];
  short Frame;
  short Remainder;
  char Strand;
  double PartialScore;
  double SRScore;
  double Score;
  pexonGFF PreviousExon;
  double GeneScore;
  char Group[MAXSTRING];
  int offset1;
  int offset2;
  short lValue;
  short rValue;
  short evidence;
  short selected;
} exonGFF;

typedef struct s_packExons             
{
  exonGFF* InitialExons;               
  exonGFF* InternalExons;
  exonGFF* TerminalExons;
  exonGFF* Singles;
  exonGFF* ORFs;

  long nInitialExons;                  
  long nInternalExons;
  long nTerminalExons;
  long nSingles;
  long nORFs;

  long nExons;                         
} packExons;

typedef struct s_packGenes
{
  exonGFF* **Ga;
  exonGFF* Ghost;
  exonGFF* GOptim;

  exonGFF* **d;
  long* km;
  long* je;
} packGenes;

typedef struct s_packEvidence
{
  site* vSites;
  long nvSites;

  exonGFF* vExons; 
  long nvExons;  

  long i1vExons;
  long i2vExons;

  long ivExons;
} packEvidence;

typedef struct s_SR
{
  long Pos1;
  long Pos2;
  double Score;
} SR;

typedef struct s_packSR
{
  SR** sRegions;
  int* nRegions;
  int nTotalRegions;
} packSR;


typedef struct s_dumpNode *pdumpNode;  
typedef struct s_dumpNode
{
   long acceptor;                      
   long donor;
   short frame;
   char strand;
   char type[MAXTYPE];

   exonGFF* exon;                      
   pdumpNode next;
} dumpNode; 

typedef struct s_dumpHash              
{
  pdumpNode* T;
  long total; 
} dumpHash;

typedef struct s_packDump
{ 
  site* dumpSites;                     
  long ndumpSites;
  
  exonGFF* dumpExons;                  
  long ndumpExons;

  dumpHash* h; 
} packDump;

typedef struct s_account               
{
  long starts, starts_r,
    stops, stops_r,
    acc, acc_r,
    don, don_r;

  long first, first_r,
    internal, internal_r,
    terminal, terminal_r,
    single, single_r,
    orf, orf_r;

  long totalExons;
  

  int tSites, tExons, tGenes, tSort, tScore, tBackup;
  time_t tStart;

} account;

typedef struct s_packGC
{
  long* GC;
  long* N;
} packGC;


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

typedef struct s_paramexons            
{
float OligoWeight;
float OligoCutoff;
float ExonWeight;
float ExonCutoff;
} paramexons;

typedef struct s_gparam                
{
  int leftValue;
  int rightValue;

  profile* StartProfile;               
  profile* AcceptorProfile;
  profile* DonorProfile;
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
(C). IMPORTED HEADERS
*************************************************************************/

void printError(char *s);

void printMess(char* s);

void printRes(char* s);

void printReadingInfo(char* s);

long GetSitesWithProfile(char *s, profile *p, site *st, long l1, long l2); 

long GetStopCodons(char *s, profile *p, site *sc, long l1, long l2);

long BuildInitialExons(site *Start, long nStarts, 
                       site *Donor, long nDonors,
                       site *Stop, long nStops,
                       int MaxDonors,
					   char* Sequence,
                       exonGFF *Exon);

long BuildInternalExons(site *Acceptor, long nAcceptors, 
                        site *Donor, long nDonors,
                        site *Stop, long nStops,
                        int MaxDonors,
						char* Sequence,
                        exonGFF* Exon);
 
long BuildTerminalExons (site *Acceptor, long nAcceptors, 
                         site *Stop, long nStops,
                         long LengthSequence,
                         long cutPoint,
						 char* Sequence,
						 exonGFF* Exon);

long BuildSingles(site *Start, long nStarts, 
                  site *Stop, long nStops,
                  long cutPoint,
				  char* Sequence,
                  exonGFF *Exon);

long BuildORFs(site *Start, long nStarts, 
			   site *Stop, long nStops,
			   long cutPoint,
			   char* Sequence,
			   exonGFF *Exon);

packSites* RequestMemorySites();
packExons* RequestMemoryExons();
exonGFF* RequestMemorySortExons();
packEvidence* RequestMemoryEvidence();
gparam* RequestMemoryParams();
packGenes* RequestMemoryGenes();
packDump* RequestMemoryDumpster();
dict* RequestMemoryAaDictionary();
void RequestMemoryProfile(profile* p);
account* RequestMemoryAccounting();

void readargv (int argc,char *argv[],
	       char *ParamFile, char* SequenceFile,
	       char *ExonsFile, char* SRFile);

int readparam (char *name, gparam** isochores);

account* InitAcc();

void OutputHeader(char* locus, long l);

int IniReadSequence(FILE* seqfile, char* line);

int ReadSequence (FILE* seqfile, char* Sequence, char* nextLocus);

long FetchSequence(char *s, char* r);

long ReadExonsGFF (char *FileName, packEvidence* pv, dict* d);

void SwitchPositions(packExons *allExons);

void SearchEvidenceExons(packEvidence* pv, long l2);

void SortExons(packExons* allExons, 
               packExons* allExons_r, 
               packEvidence* pv,
               exonGFF* Exons,         
               long l1, long l2);

void SwitchCounters(packEvidence* pv);

void Output(packSites* allSites, packSites* allSites_r,
            packExons* allExons, packExons* allExons_r,
            exonGFF* exons, long nExons, char* Locus, 
	    long l1, long l2, char* Sequence, gparam* gp, dict* dAA);

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
                char* Sequence, gparam* gp, dict* dAA);

void OutputStats(char* Locus);
void OutputTime();

void RecomputePositions(packSites* allSites, long l);

void cleanAcc(account* m);

void PrintProfile (profile *p, char* signal);

long ReadGeneModel (FILE *file, dict *d, int nc[], int ne[], 
                    int UC[][MAXENTRY], int DE[][MAXENTRY], 
		    long md[], long Md[], int block[]);

void PrintSites (site *s, long ns,int type,
                 char Name[], int Strand,
                 long l1, long l2,
                 char* seq,
                 profile *p);

void PrintExons (exonGFF *e, long ne, int type, char Name[],
                 long l1, long l2, char* Sequence,  dict* dAA); 

void resetDict(dict* d);

int setkeyDict(dict *d, char s[]);

int getkeyDict(dict *d, char s[]);

int Translate(long p1, long p2, short fra, short rmd,
              char *s, dict* dAA, char sAux[]);

void ReverseSubSequence(long p1, long p2, char* s, char* r);

void CorrectExon(exonGFF *e);

void CorrectORF(exonGFF *e);

void SwitchFrames(exonGFF* e, long n);

void SwitchFramesDa(packGenes* pg, int nclass);

void SwitchFramesDb(packGenes* pg, int nclass);

void UndoFrames(exonGFF *e, long n);

void BuildSort(dict *D, int nc[], int ne[], int UC[][MAXENTRY], 
               int DE[][MAXENTRY], int nclass, long km[], 
	       exonGFF* **d, exonGFF *E, long nexons);

void PrintSite(site *s, int type, char Name[], int Strand,
               char* seq, profile *p);

void PrintGExon(exonGFF *e, char Name[], char* s, dict* dAA, 
               long ngen, int AA1, int AA2, int nAA);

void PrintXMLExon(exonGFF *e, char Name[], 
		  long ngen, long nExon, 
		  int type1, int type2, 
		  int nExons);

void TranslateGene(exonGFF* e,
                   char* s,
                   dict* dAA,
                   long nExons,
                   int** tAA,
                   char* prot,
                   long* nAA);

void resetDumpHash(dumpHash* h);

void setExonDumpHash(exonGFF* E, dumpHash *h);

exonGFF* getExonDumpHash(exonGFF* E, dumpHash *h);
 
void setAADict(dict *d, char s[], char aA);

char getAADict(dict *d, char s[]);

void CookingGenes(exonGFF *e, char Name[], char* s,
                  gparam* gp, dict* dAA);

float MeasureSequence(long l1,long l2,char* s);

packSR* RequestMemorySimilarityRegions();

gparam ** RequestMemoryIsochoresParams();

long ReadSR (char *FileName, packSR* sr, long LengthSequence);

void shareGeneModel(gparam** isochores, int n);

long OligoToInt(char *s, int ls, int cardinal);

char* RequestMemorySequence(long L);

void ScoreExons(char *Sequence, 
                packExons* allExons, 
                long l1,
                long l2,
                int Strand,
                packSR* sr,
                gparam** isochores,
                int nIsochores,
                packGC* GCInfo);

void GetcDNA(exonGFF* e, char* s, long nExons, char* cDNA, long* nNN);

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
			  long l1, long l2,
			  int Strand,
			  packSR* sr,
			  gparam* gp,
			  gparam** isochores,
			  int nIsochores,
			  packGC* GCInfo);

void resetEvidenceCounters(packEvidence* pv);

void ComputeStopInfo(exonGFF* e, char* s);
