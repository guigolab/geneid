/***************************************************************
*                                                              *
*   BLAST2GFF.                                                 *
*   (c) by Enrique Blanco &                                    *
*   Roderic Guigo, 2000                                        *
*                                                              *
*   Module: blast2gff.h                                        *
*                                                              *
*   Definition of all data structures used by blast2gff        *
***************************************************************/

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

/* Definitions used in blast2gff */
#define NUMHSPS 6000000                /* maximum number of HSPs          */
#define NUMSR 5000000                  /* maximum number of SRs           */
/* Those values were increased for chromosome analysis */

#define LOCUSLENGTH 50                 /* maximum number of chars (locus) */

#define VERSION   "blast2gff"          /* The name of the game            */
#define HSP       "HSP"                /* The name of the sites           */
#define NONAME    "NONAME"
#define SR        "SR"                

#define FRAMES 3                       /* 3 possible reading frames       */
#define STRANDS 2                      /* 2 strands of DNA molecule       */ 

#define FORWARD 0                      /* DNA - Strands                   */
#define REVERSE 1

#define NOTFOUND -1

#define COFFSET 1                      /* C-arrays correction             */
#define LENGTHCODON 3                  /* 3 nucleotids                    */

#define MAXLINE 1000                   /* Max number of chars/inputline   */

#define MAXTYPE 20                     /* Max number of chars/exonTypes   */

#define MAXSTRING 100                   

#define FILENAMELENGTH 200             /* maximum length of filenames     */
#define INFI 9999999                   /* the biggest number of the world */
#define MAXSCORE 10000.0               /* Infinity score (evidence exons) */

#define MIN(a,b) (a<b)?a:b;
#define MAX(a,b) (a>b)?a:b;

#define MINUTE 60                      /* 1 minute = 60 seconds           */     

/****************************************************************************/

/* GeneId data types */
typedef struct s_hsp               /* Blast HSP           */
{
  long StartQ;
  long EndQ;
  char StrandQ;
  short FrameQ;
  float Score;
  float Bits;
  double Probability;
  double Expect;
  int Similarity;
  int Identity;
  long StartS;
  long EndS;
  long Lenght;
  char Subject[LOCUSLENGTH];
} hsp;

typedef struct s_packHSP             /* Pack of hsps                   */
{
  hsp** *hsps;
  long* nHsps;
  long nTotalHsps;
} packHSP;     

typedef struct s_sr_t* psr_t;
typedef struct s_sr_t 
{
  long Start;
  long End;
  char Strand;
  short Frame;
  float Score;
  char Locus[LOCUSLENGTH];
  psr_t next;
} sr_t;

typedef struct s_packSR
{
  sr_t* sr[STRANDS*FRAMES];            /* Lists of SR */
  int nr[STRANDS*FRAMES];              /* Counter of every list */
  sr_t* nr0[STRANDS*FRAMES];           /* Pointer to processed part of every list */

  sr_t* sortSR[NUMSR];
  long nSR;
} packSR;

typedef struct s_account               /* Accounting                      */
{
  time_t tStart;
} account;                      

/*** Headers of functions ***/
void printMess(char* s);
void printError(char *s);
void OutputTime();
void readargv (int argc,char* argv[],char* HSPFile);
void ReadHSP (packHSP* allHsp,char* HSPFile, char Query[LOCUSLENGTH]);
packHSP* RequestMemoryHSP();
packSR* RequestMemorySR();
account* InitAcc();
void Output(packHSP* allHSP, packSR* allSR, char Query[LOCUSLENGTH]);
void ProjectHSP(packHSP *allHsp, packSR *allSr);
void SortHSP(packHSP *allHsp);
void JoinSR(packSR* allSr);
void SortSR(packSR *allSr);
