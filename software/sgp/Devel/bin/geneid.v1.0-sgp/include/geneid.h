/*************************************************************************
*                                                                        *
*   Module: geneid.h                                                     *
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

/* $Id: geneid.h,v 1.2 2000-07-21 11:40:12 eblanco Exp $ */

/* Include libraries */
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

/* Definitions used in GeneId */

#define NUMSITES 25000                 /* maximum number of sites         */
#define NUMEXONS 50000                 /* maximum number of exons         */

#define NUMEEVIDENCES 10000            /* maximum number of evidence exons */
#define NUMSEVIDENCES 2*NUMEEVIDENCES          

#define MAXGENE 10000                  /* Max number of genes(multigenes) */
#define MAXEXONGEN 50000               /* Max numer of exon per gene      */
#define LENGTHSi 100000                /* Max length of each split of S   */ 
#define OVERLAP 15000                  /* Max length of a exon(overlap)   */

#define MAXBACKUPSITES 125000          /* Backup Information sizes */
#define MAXBACKUPEXONS 125000      
#define MAXDUMPHASH 124997

#define MAXAA 10000                    /* Max lenght (aminoacids/protein) */

#define LOCUSLENGTH 50                 /* maximum number of chars (locus) */

#define PROFILEDIM 100                 /* maximum profile's dimension     */

#define OLIGOLENGTH 10                 /* maximum oligo-word length       */

#define EXONLENGTH 18                  /* minimum exon length (internal   */
                                       
#define SINGLEGENELENGTH 60            /* minimum single gene length      */ 

#define INTRONLENGTH 50                /* minimum exon length (internal   */
                                       /* and terminal)                   */
#define NULL_OLIGO_SCORE  -4           /* temporal hack to recompute */

#define VERSION   "geneid_v1.0"        /* The name of the game            */
#define SITES     "geneid_v1.0"        /* The name of the sites           */
#define EXONS     "geneid_v1.0"        /* The name of the exons           */
#define EVIDENCE  "evidence"           /* The name of evidence exons      */

#define FRAMES 3                       /* 3 possible reading frames       */
#define STRANDS 2                      /* 2 strands of DNA molecule       */ 

#define FORWARD 0                      /* DNA - Strands                   */
#define REVERSE 1

#define ACC 0                          /* Define site types               */
#define DON 1
#define STA 2
#define STO 3

#define FIRST    0                     /* Define exon types               */
#define INTERNAL 1
#define TERMINAL 2
#define SINGLE   3

#define sACC "Acceptor"                /* Define site types               */
#define sDON "Donor"
#define sSTA "Start"
#define sSTO "Stop"

#define sFIRST    "First"              /* Define exon types               */
#define sINTERNAL "Internal"
#define sTERMINAL "Terminal"
#define sSINGLE   "Single"

#define BLOCK 1                        /* define a block rule             */
#define NONBLOCK 0                     /* define a non/block rule         */

#define COFFSET 1                      /* C-arrays correction             */
#define LENGTHCODON 3                  /* 3 nucleotids                    */

#define MAXLINE 1000                   /* Max number of chars/inputline   */

#define MAXTYPE 20                     /* Max number of chars/exonTypes   */

#define MAXENTRY 97                    /* Dictionary definitions          */
#define MAXSTRING 100                   
#define MAXINFO 100
#define NOTFOUND -1

#define MEGABYTE 1048576               /* Accounting                      */
#define MAXTIMES 100

#define PARAMETERFILE  "param.default" /* default parameters file         */
					

#define FILENAMELENGTH 200             /* maximum length of filenames     */
#define INF 1.7976931348623157E+308    /* maximal decimal value of double */
                                       /* from limits.h                   */
#define INFI 9999999                   /* the biggest number of the world */
#define MAXSCORE 10000.0               /* Infinity score (evidence exons) */

#define MINUTE 60                      /* 1 minute = 60 seconds           */

#define MIN(a,b) (a<b)?a:b;
#define MAX(a,b) (a>b)?a:b;

#define MAXSR 1000                     /* Maximum number of sr per frame  */
#define MAXISOCHORES 4                 /* Maximum number of isochores */
#define PERCENT 100   

#define NOVALUE 0
#define NOGROUP -13

#define ISOCONTEXT 1000

/****************************************************************************/

/* GeneId data types */
typedef struct s_profile               /* Predict sites profile           */
{
  int    dimension;
  int    offset;
  float  cutoff;
  int    order;

  long dimensionTrans;
  float*  transitionValues[PROFILEDIM];
  /*  float  matrix[127][PROFILEDIM]; */
} profile;

typedef struct s_site                  /* Site data type                  */
{
  long Position;                       
  double Score;                        
} site;

typedef struct s_packSites             /* Pack of sites                   */
{
  site* StartCodons;                   /* Sites                           */
  site* AcceptorSites;
  site* DonorSites;
  site* StopCodons;

  long  nStartCodons;                  /* Counters of sites               */
  long  nAcceptorSites;
  long  nDonorSites;
  long  nStopCodons;

  long nSites;                         /* Counter of the total number     */
} packSites;

typedef struct s_exonGFF *pexonGFF;    /* Standard GFF-format exon        */
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
  int Group;
  int offset1;
  int offset2;
  short evidence;
  short selected;
} exonGFF;

typedef struct s_packExons             /* Pack of exons                   */
{
  exonGFF* InitialExons;               /* Exons                           */
  exonGFF* InternalExons;
  exonGFF* TerminalExons;
  exonGFF* Singles;

  long nInitialExons;                  /* Counters of exons               */
  long nInternalExons;
  long nTerminalExons;
  long nSingles;

  long nExons;                         /* Counter of the total number     */
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
  int* iRegions;
  int* nRegions;
  int nTotalRegions;
} packSR;


typedef struct s_dumpNode *pdumpNode;  /* Structure of a sinonimous       */
typedef struct s_dumpNode
{
   long acceptor;                      /* Hash keys                       */
   long donor;
   short frame;
   char strand;
   char type[MAXTYPE];

   exonGFF* exon;                      /* Info                            */
   pdumpNode next;
} dumpNode; 

typedef struct s_dumpHash              /* Structure of a hash table       */
{
  pdumpNode T[MAXDUMPHASH];
  long total; 
} dumpHash;

typedef struct s_packDump
{ 
  site* dumpSites;                     /* backup sites                    */
  long ndumpSites;
  
  exonGFF* dumpExons;                  /* backup exons                    */
  long ndumpExons;

  dumpHash* h; 
} packDump;

typedef struct s_account               /* Accounting                      */
{
  long starts, starts_r,
       stops, stops_r,
       acc, acc_r,
       don, don_r;

  long first, first_r,
       internal, internal_r,
       terminal, terminal_r,
       single, single_r;

  long totalExons;
  

  int tSites, tExons, tGenes, tSort, tScore, tBackup;
  time_t tStart;

} account;

typedef struct s_node *pnode;          /* Structure of a sinonimous       */
typedef struct s_node
{
  char s[MAXSTRING];
  int key;
  pnode next;
} node; 

typedef struct s_dict                  /* Structure of a hash table       */
{
  pnode T[MAXENTRY];
  int nextFree; 
} dict;

typedef struct s_paramexons            /* Data type of score-params       */
{
float OligoWeight;
float OligoCutoff;
float ExonWeight;
float ExonCutoff;
} paramexons;

typedef struct s_gparam                /* Data type of GeneId-params      */
{
  int leftValue;
  int rightValue;

  profile* StartProfile;               /* Profiles for predicting sites   */
  profile* AcceptorProfile;
  profile* DonorProfile;
  profile* StopProfile;

  float* OligoLogsIni[3];              /* Markov initial values           */
  float* OligoLogsTran[3];

  long OligoDim;                       /* Markov models parameters        */
  long OligoDim_1;           
  int  OligoLength;

  paramexons* Initial;                 /* Score exons parameters          */
  paramexons* Internal;
  paramexons* Terminal;
  paramexons* Single;
                                                                            
  int MaxDonors;                       /* Max number of exons per acc     */

  dict* D;                             /* GeneModel configuration         */
  int*  nc;
  int*  ne;
  long* md;
  long* Md;
  int UC[MAXENTRY][MAXENTRY];
  int DE[MAXENTRY][MAXENTRY];
  int block[MAXENTRY];
  int nclass;

} gparam;


/* Import function headers */

/* Printing error messages */
void printError(char *s);

/* Printing messages (information) */
void printMess(char* s);

/* Printing partial and final results (information) */
void printRes(char* s);

long GetSitesWithProfile(char *s, profile *p, site *st, long l1, long l2); 

long GetStopCodons(char *s, profile *p, site *sc, long l1, long l2);

long BuildInitialExons(site *Start, long nStarts, 
                       site *Donor, long nDonors,
                       site *Stop, long nStops,
                       int MaxDonors,
                       exonGFF *Exon);

long BuildInternalExons(site *Acceptor, long nAcceptors, 
                        site *Donor, long nDonors,
                        site *Stop, long nStops,
                        int MaxDonors,
                        exonGFF* Exon);
 
long BuildTerminalExons (site *Acceptor, long nAcceptors, 
                         site *Stop, long nStops,
                         long LengthSequence,
                         exonGFF* Exon,
                         long cutPoint);

long BuildSingles(site *Start, long nStarts, 
                  site *Stop, long nStops,
                  long cutPoint,
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

long SearchEvidenceExons(packEvidence* pv, long l2);

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

void OutputStats(char* Locus, long nvExons);
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

void SwitchFrames(exonGFF *e, long n, int when);

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

void TranslateGen(exonGFF* e, char* s,dict* dAA, long nExons, 
                  int tAA[MAXEXONGEN][2], char* prot, int* nAA);

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

long OligoToInt(char *s, int ls);

char* RequestMemorySequence(long L);


void ScoreExons(char *Sequence, 
		packExons* allExons, 
		long l1, long l2, int Strand, packSR* sr,
		gparam** isochores, int nIsochores);
