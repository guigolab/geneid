package SGP::filters;
use strict;
#
#line 462 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
######################################################################
#                         SGPfilters.pm                              #
######################################################################
#
#     Modules to filter HSPs and project them to SRs.
#
#     Copyright (C) 2001 - Josep Francesc ABRIL FERRANDO
#                                   Genis PARRA FARRE
#                                 Enrique BLANCO GARCIA
#                                 Roderic GUIGO SERRA
#line 1653 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
######################################################################
#
#line 473 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
# 
#line 1677 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
# $Id: SGPfilters.pm,v 1.1 2001-06-05 15:25:35 jabril Exp $
#line 430 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
#
require Exporter;
@SGP::filters::ISA = qw(Exporter);
@SGP::filters::EXPORT = qw(
                            
                            
                            );
$SGP::filters::VERSION = '0.1';

use Inline ( C => DATA,
             NAME => 'SGP::filters',
             VERSION => '0.10' );

__DATA__

__C__

#line 506 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
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
#define NUMHSPS 50000                  /* maximum number of HSPs          */
#define NUMSR 50000                    /* maximum number of SRs           */

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

#line 638 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
/* Setup Flags of blast2gff */
int   VRB=0, GFFIN, PSR=1 ,GFFOUT=0;

/* Accounting time and results */
account *m;  

int main (int argc, char *argv[])
{
  /* Input Files */
  char   HSPFile[FILENAMELENGTH];
  
  /* HSP data structures */
  packHSP* allHsp;

  /* SR data structures */
  packSR* allSr;

  /* Query name */
  char Query[LOCUSLENGTH];

  char mess[MAXSTRING];

  /* initializing stats... */
  m = (account*)InitAcc(); 

  /* Read setup options */
  readargv(argc,argv,HSPFile);

  /* Alloc main memory structures */
  printMess("Request Memory to System");
  allHsp =      (packHSP*) RequestMemoryHSP();
  allSr  =      (packSR*) RequestMemorySR(); 

  /* Reading HSPs */
  printMess("Reading HSP file");
  ReadHSP(allHsp,HSPFile,Query);
  sprintf(mess,"%ld HSPs read",allHsp->nTotalHsps);
  printMess(mess);
  
  /* Sorting HSPs */
  printMess("QuickSorting HSPs");
  SortHSP(allHsp);
  
  /* Project HSP to SR */
  if (PSR)
    {
      printMess("Project HSPs to SRs");
      ProjectHSP(allHsp, allSr);

      printMess("Simplifiying SRs");
      JoinSR(allSr);

      printMess("MergeSorting SRs");
      SortSR(allSr);
    }
  
  /* Printing Results */
  printMess("Printing selected data");
  Output(allHsp, allSr, Query);

  /* The End */
  OutputTime(); 
  exit(0);
  return(0);   
}

#line 707 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
extern int VRB, GFFIN, GFFOUT, PSR;
                
extern char* optarg;
extern int optind;

char *USAGE="Incorrect usage:\nNAME\n\tblast2gff - A program to build similarity regions from blastx outputs\nSYNOPSIS\n\tblast2gff <HSP_file>\n\n";

void printHelp()
{
  printf("Short Manual for blast2gff:\n");
  printf("------------------------\n\n");
  printf("Setup Options:\n\n");
 
  printf("\t-g: Input HSP in GFF format\n");
  printf("\t-G: Print HSP in GFF format\n");
  printf("\t-o: Not SR construction\n");
  printf("\t-v: Verbose\n");
  printf("\t-h: Show this Short Manual\n");
}

void readargv (int argc,char* argv[],
	       char* HSPFile) 
{
  int c;
  int error=0;
  char mess[MAXSTRING];

  /* Reading setup options */
  while ((c = getopt(argc,argv,"Ggvoh")) != -1)
    switch(c)
      { 
      case 'g': GFFIN++; 
	break;
      case 'G': GFFOUT++; 
	break;
      case 'o': PSR--; 
	break;
      case 'v': VRB++; 
	break;
      case 'h': printHelp();
	exit(1);
	break;
      }

  if (error)
    printError(USAGE);

  /* Setup Errors: Wrong number of filenames */
  /* Get the name of the file *.hsp */
  if (optind < argc)
    {
      strcpy(HSPFile,argv[optind]);
      optind++;
      if (optind < argc)
	{ 
	  sprintf(mess,"Too many files. Only one filename needed\n%s",USAGE);
	  printError(mess);
	}
    }
  else
    {
      sprintf(mess,"Where is the HSP file?\n%s",USAGE);
      printError(mess);
    }
}

#line 776 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
/* Init accounting */
account* InitAcc()
{
  account* m; 
  
  m = (account *) malloc(sizeof(account));

  /* Reset counters */

  printMess("Reset accounting");
  /* What time is it? */ 
  clock();
  (void) time(&m->tStart);

  return(m);
}

#line 796 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
packHSP* RequestMemoryHSP()
{
  packHSP *allHsp;
  int i; 

  if ((allHsp =
       (packHSP*) malloc(sizeof(packHSP)))  == NULL)
    printError("Not enough space to hold HSP data structure");

  /* HSP */
  if ((allHsp->hsps =
       (hsp ** *) calloc(STRANDS*FRAMES, sizeof(hsp**))) == NULL)
    printError("Not enough space to hold HSP 6-array");  


  for(i=0;i<STRANDS*FRAMES;i++)
    if ((allHsp->hsps[i] =
	 (hsp**) calloc(NUMHSPS, sizeof(hsp*)))  == NULL)
      printError("Not enough space to hold HSPs");

  /* Counters */
  if ((allHsp->nHsps =
       (long*) calloc(STRANDS*FRAMES, sizeof(int)))  == NULL)
    printError("Not enough space to hold HSP counters");  

  /* Hack of Alpha */
  for(i=0;i<STRANDS*FRAMES;i++)
    allHsp->nHsps[i] = 0;

  allHsp->nTotalHsps = 0;

  return(allHsp);
}

packSR* RequestMemorySR()
{
  packSR *allSr;
  int i;

  if ((allSr =
       (packSR*) malloc(sizeof(packSR)))  == NULL)
    printError("Not enough space to hold SR data structure");

  for(i=0;i<STRANDS*FRAMES;i++)
    {
      allSr->nr[i]  = 0;
      allSr->nr0[i] = NULL;
      allSr->sr[i]  = NULL;
    }
  
  allSr->nSR = 0;
  
  return(allSr);
}

#line 854 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
extern int GFFIN;

hsp* allocateNewHSP()
{
  hsp* p;
  
  /* New HSP */
  if ((p = (hsp *) malloc(sizeof(hsp))) == NULL)
    printError("Not enough space to hold one new HSP");
  
  return(p);
}

long ReadHSP_RAW(packHSP* h,char *FileName)
{
  int i;
  FILE *file;
  char line[MAXLINE];
  char mess[MAXSTRING];
  long pos1, pos2;
  int score;
  char strand;
  short frame;
  float bits;
  double prob, exp;
  int id,sim;
  long pos3,pos4;
  int indexFrame;
  long L;
  char Subject[LOCUSLENGTH];

  if ((file=fopen(FileName, "r"))==NULL)
    printError("The HSP file file cannot be open for read");
  
  i = 0;
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      if(line[0]=='#')
	{
	  /* Skip this line */
	}
      else
	{
	  if ((sscanf(line, "%*s %*s %ld %ld %c %hd %d %f %lf %lf %d %d %s %ld %ld %*c %*d %ld",
		      &pos1,
		      &pos2,
		      &strand,
		      &frame,
		      &score,
		      &bits,
		      &prob,
		      &exp,
		      &id,
		      &sim,
		      Subject,
		      &pos3,
		      &pos4,
		      &L) != 14) || (frame<1) || (frame>3))
	    {
	      sprintf(mess, "Error reading HPS: line %d\n",i);
	      printError(mess);
	    }
	  else
	    {
	      /* Allocating HSPs in packHSP */
	      if (strand == '+')
		indexFrame = (frame % 3);
	      else
		indexFrame = (frame % 3) + FRAMES;

	      /* New item HSP */
	      h->hsps[indexFrame][h->nHsps[indexFrame]] = allocateNewHSP();

	      /* Setting attributes */
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StartQ = pos1;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->EndQ = pos2;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StrandQ = strand;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->FrameQ = frame;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Score = score;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Bits = bits;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Probability = prob;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Expect = exp;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Similarity = sim;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Identity = id;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StartS = pos3;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->EndS = pos4;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Lenght = L;
	      strcpy(h->hsps[indexFrame][h->nHsps[indexFrame]]->Subject,Subject);
	      h->nHsps[indexFrame]++;
	      h->nTotalHsps++;
	    }
	}
      i++;
    }
  fclose(file);
  
  return(h->nTotalHsps);
}

long ReadHSP_GFF (packHSP* h,char *FileName, char Query[LOCUSLENGTH])
{
  int i;
  FILE *file;
  char line[MAXLINE];
  char mess[MAXSTRING];
  long pos1, pos2;
  char strand;
  short frame;
  float score;
  int indexFrame;
  char Subject[LOCUSLENGTH];


  if ((file=fopen(FileName, "r"))==NULL)
    printError("The HSP file file cannot be open for read");
  
  i = 0;
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      if(line[0]=='#')
	{
	  /* Skip this line */
	}
      else
	{
	  if ((sscanf(line, "%s %*s %*s %ld %ld %f %c %hd %s",
		      Query,
		      &pos1,
		      &pos2,
		      &score,
		      &strand,
		      &frame,
		      Subject) != 7) || (frame<1) || (frame>3))
	    {
	      sprintf(mess, "Error reading HPS: line %d\n",i);
	      printError(mess);
	    }
	  else
	    {
	      /* Allocating HSPs in packHSP */
	      if (strand == '+')
		indexFrame = (frame % 3);
	      else
		indexFrame = (frame % 3) + FRAMES;

	      /* New item HSP */
	      h->hsps[indexFrame][h->nHsps[indexFrame]] = allocateNewHSP();
	      
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StartQ = pos1;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->EndQ = pos2;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->StrandQ = strand;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->FrameQ = frame;
	      h->hsps[indexFrame][h->nHsps[indexFrame]]->Score = score;
	      strcpy(h->hsps[indexFrame][h->nHsps[indexFrame]]->Subject,Subject);

	      h->nHsps[indexFrame]++;
	      h->nTotalHsps++;
	    }
	}
      i++;
    }
  fclose(file);
  
  return(h->nTotalHsps);
}

void ReadHSP (packHSP* allHsp,char* HSPFile, char Query[LOCUSLENGTH])
{
  if (GFFIN)
    ReadHSP_GFF(allHsp,HSPFile,Query); 
  else
    {
      ReadHSP_RAW(allHsp,HSPFile); 
      strcpy(Query,NONAME);
    }
}

#line 1034 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
int Split(hsp** hsps,int i, int j)
{
  long k;
  int x,c1,c2,m;
  int pivot;
  int pivotStartQ;
  hsp* tmp;
  
  /* Choosing the random pivot */
  k = hsps[i]->StartQ;

  /* How many elements have less StartQ value than pivot */
  for(m = 0, x=i+1; x <=j; x++)
    if (hsps[x]->StartQ <= k)
	m++;
  
  /* This is the right place for the pivot */
  pivot = m+i;

  tmp = hsps[i];
  hsps[i] = hsps[pivot];
  hsps[pivot] = tmp;

  pivotStartQ = hsps[pivot]->StartQ; 
  c1 = i;
  c2 = pivot+1;
  while (c1 < m+i && c2 <= j) 
    if (hsps[c1]->StartQ <= pivotStartQ)
	c1++;
    else
      {
	tmp = hsps[c1];
	hsps[c1] = hsps[c2];
	hsps[c2] = tmp;

	c2++;
      }

  return(pivot);
}

void quickSort(hsp** hsps, int i, int j)
{
  int pivot;

  if (i == j+1)
    /* Nothing */;
  else
    {
      pivot = Split(hsps,i,j);
      
      quickSort(hsps,i,pivot-1);
      quickSort(hsps,pivot+1,j);
    }
}

void SortHSP(packHSP *allHsp)
{
  int i;

  /* quicksorting all the strands and frames */ 
  for(i=0; i<STRANDS*FRAMES; i++)
    quickSort(allHsp->hsps[i],0,allHsp->nHsps[i]-1);
}

#line 1102 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
void printSR(sr_t* t, long i)
{
  sr_t* tmp; 

  printf("Resultados iteracion %ld\n",i);
  
  tmp = t;
  while(tmp != NULL)
    {
      printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\n",
	     NONAME,
	     VERSION,
	     SR,
	     tmp->Start,
	     tmp->End,
	     tmp->Score,
	     tmp->Strand,
	     tmp->Frame); 
      tmp = tmp->next;
    }
}

int FreeItems(sr_t* q)
{
  int numElems;
  
  if (q == NULL)
    numElems = 0;
  else
    {
      numElems = FreeItems(q->next);
      numElems++;
      free(q);
    }
  return(numElems);
}

/* Returns the number of free SR */
int CutAndPaste(sr_t** Left, sr_t** Right,
                sr_t** insLeft, sr_t** insRight,
                sr_t** sr)
{
  sr_t* tmp;
  int numElems;
  
/*    printf("nr0: %ld %ld %d\n",Left->Start,Left->End,Left->Score); */
/*    printf("sr: %ld %ld %d\n",sr->Start,sr->End,sr->Score); */
/*    printf("firstSR: %ld %ld %d\n",Right->Start,Right->End,Right->Score); */
/*    printf("fnSR: %ld %ld %d\n",insLeft->Start,insLeft->End,insLeft->Score); */
/*    printf("newSR: %ld %ld %d\n",insRight->Start,insRight->End,insRight->Score); */
  
  /* Cut and Paste at the beginning of the list */
  if (*Left == NULL)
    {
      /* tmp used for freeing the useless nodes*/
      tmp = *sr;
      *sr = *insLeft;
    }
  else
    {
      /* tmp used for freeing the useless nodes*/
      tmp = (*Left)->next;
      (*Left)->next = *insLeft;
    }

  
  (*insRight)->next = (*Right)->next;
  
  /* Free nodes between [Left->next | sr] and Right */
  (*Right)->next = NULL;

  numElems = FreeItems(tmp);
  
  return(numElems);
}

void Chain(sr_t** newSR, sr_t* auxSR, sr_t** firstNewSR)
{
  /* Chaining the new SR to the chain */
  if (*newSR != NULL)
    {
      (*newSR)->next = auxSR;
      (*newSR) = (*newSR)->next;
    }
  else
    {
      /*        printMess("INT: iniciando la nueva cadena"); */
      
      *newSR = auxSR;
      *firstNewSR = auxSR;
    }
}

sr_t* RequestMemoryNewSR()
{
  sr_t* s;
  
  if ((s = (sr_t *) malloc(sizeof(sr_t))) == NULL)
    printError("Not enough space to hold one new SR");
  
  return(s);
}

void project1HSP(hsp** hsps, long nHSP, packSR *allSr)
{
  int index;
  long i;
  sr_t* beforeFirstSR;
  sr_t* firstSR;
  sr_t* auxSR;
  sr_t* newSR;
  sr_t* firstNewSR;
  int numElems;
  
  /* For each HSP do... */
  for(i=0; i<nHSP; i++)
    {
      /*    printf("\nTratando HSP %ld: %ld %ld %d %hd\n", */
      /*  	     i,hsps[i]->StartQ,hsps[i]->EndQ,hsps[i]->Score,hsps[i]->FrameQ); */
      
      index = (hsps[i]->StrandQ == '+')?
      hsps[i]->FrameQ % 3 :  
	(hsps[i]->FrameQ % 3) + FRAMES;
      
      /* Searching into the list, the first SR after this HP */
      /* beforeFirstSR: the last before the first SR after this HP */
      /* Possible values for firstSR and beforeFirstSR */
      /* a. Reading from the start (1 element) */
      if (allSr->nr[index] > 0 && allSr->nr0[index]==NULL)
	{
	  firstSR = allSr->sr[index];
	  beforeFirstSR = NULL;
	}
      else
	{
	  /* b. Empty list */
	  if (allSr->nr[index] == 0 && allSr->nr0[index] == NULL)         
	    firstSR = NULL;         
	  /* c. Normal situation: inside the list */
	  else
	    firstSR = allSr->nr0[index];
	  
	  beforeFirstSR = firstSR;
	}
      
      /* scanning the list of SR */
      while(firstSR != NULL && firstSR->End < hsps[i]->StartQ)
	{
	  beforeFirstSR = firstSR;
	  firstSR = firstSR->next;
	}
      
      /* last SR before first selected SR */
      allSr->nr0[index] = beforeFirstSR;
      
      /* Empty list: Insert first SR (the same HP) */
      if(firstSR == NULL && beforeFirstSR == NULL)
	{
	  /*  printMess("Lista vacia: creando primer SR\n"); */
	  
	  /* Creating a new SR placed after end of this list */
	  auxSR = RequestMemoryNewSR();
	  
	  /* Setting HSP values */
	  auxSR->Start = hsps[i]->StartQ;
	  auxSR->End = hsps[i]->EndQ;
	  auxSR->Strand = hsps[i]->StrandQ;
	  auxSR->Frame = hsps[i]->FrameQ;
	  auxSR->Score = hsps[i]->Score;
	  strcpy(auxSR->Locus,hsps[i]->Subject);
	  auxSR->next = NULL;
	  
	  /* Insert at the beginning */
	  allSr->sr[index] = auxSR;
	  allSr->nr[index]  = 1;
	}
      else
	/* Insert SR at the end of the list */
	if(firstSR == NULL && beforeFirstSR != NULL)
	  {
	    /*  printMess("Recorrida lista completa: insertar al final\n"); */
	    
	    /* Creating a new SR placed after end of this list */
	    auxSR = RequestMemoryNewSR();
	    
	    /* Setting HSP values */
	    auxSR->Start = hsps[i]->StartQ;
	    auxSR->End = hsps[i]->EndQ;
	    auxSR->Strand = hsps[i]->StrandQ;
	    auxSR->Frame = hsps[i]->FrameQ;
	    auxSR->Score = hsps[i]->Score;
	    strcpy(auxSR->Locus,hsps[i]->Subject);
	    auxSR->next = NULL;
	    
	    /* Insert at the end */
	    beforeFirstSR->next = auxSR;
	    allSr->nr[index]++;
	  }
	else
	  {
	    /* We work with the sublist starting at firstSR */
	    /* newSR is the current created node */
	    /* firstNewSR is the beginning of the new SR chain */
	    newSR = NULL;
	    firstNewSR = NULL;
	    
	    /*  printMess("Intersecciones"); */
	    /*  	    printf("INT: SR1: %ld %ld\n",firstSR->Start,firstSR->End); */
	    
	    /* 1. Creating a new SR: START(firstSR) - START(HSP) */
	    if (firstSR->Start < hsps[i]->StartQ)
	      {
		newSR = RequestMemoryNewSR();
		firstNewSR = newSR;
		
		/* Setting HSP values */
		newSR->Start = firstSR->Start;
		newSR->End = hsps[i]->StartQ-1;
		newSR->Strand = firstSR->Strand;
		newSR->Frame = firstSR->Frame;
		newSR->Score = firstSR->Score;
		strcpy(newSR->Locus,hsps[i]->Subject);
		newSR->next = NULL; 
		
		allSr->nr[index]++;
		
		/*  printf("INT: Generado primer bloque\n"); */
		/*  		printf("%ld %ld\n",newSR->Start,newSR->End); */
	      }
	    
	    /* 2. Creating the rest of intersections until last*/
	    auxSR = RequestMemoryNewSR();
	    auxSR->Start = hsps[i]->StartQ;
	    
	    while((hsps[i]->EndQ > firstSR->End) && (firstSR->next != NULL))
	      {
		/*  printf("INT: Dentro del while\n"); */
		auxSR->End = firstSR->End;
		auxSR->Score = MAX(firstSR->Score,hsps[i]->Score);
		auxSR->Strand = firstSR->Strand;
		auxSR->Frame = firstSR->Frame;
		/*  auxSR->Locus = strcat(firstSR->Locus,hsps[i]->Subject); */
		strcpy(auxSR->Locus,hsps[i]->Subject);
		auxSR->next = NULL;
		
		allSr->nr[index]++;
		
		/*  printf("INT: Creado nuevo SR(while)\n"); */
		
		Chain(&newSR,auxSR,&firstNewSR);
		
		/* jumping to the next SR on the chain */
		firstSR = firstSR->next;
		
		auxSR = RequestMemoryNewSR();
		auxSR->Start = firstSR->Start;
	      }/* endwhile */
	    
	    /* 3. Creating the last intersections */
	    /* NO Sobresale por detras el SR */
	    if (hsps[i]->EndQ >= firstSR->End)
	      {
		auxSR->End = firstSR->End;
		auxSR->Score = MAX(firstSR->Score,hsps[i]->Score);
		auxSR->Strand = firstSR->Strand;
		auxSR->Frame = firstSR->Frame;
		/*  auxSR->Locus = strcat(firstSR->Locus,hsps[i]->Subject); */
		strcpy(auxSR->Locus,hsps[i]->Subject);
		auxSR->next = NULL;
		
		/*  printf("auxSR %ld %ld %d\n",auxSR->Start,auxSR->End,auxSR->Score); */

		allSr->nr[index]++;
		
		Chain(&newSR, auxSR,&firstNewSR); 

		/* Coinciden los finales */
		if (hsps[i]->EndQ == firstSR->End)
		  {}/*  printf("INT2: igual final. No se crean mas\n"); */
		else
		  {
		    /* Creating a new SR: END(firstSR) - END(HSP) */
 		    auxSR = RequestMemoryNewSR();
		    auxSR->Start = firstSR->End+1;
		    auxSR->End = hsps[i]->EndQ;
		    auxSR->Strand = firstSR->Strand;
		    auxSR->Frame = firstSR->Frame;
		    auxSR->Score = hsps[i]->Score;
		    strcpy(auxSR->Locus,hsps[i]->Subject);
		    auxSR->next = NULL;
		    
		    allSr->nr[index]++;
		    
		    Chain(&newSR,auxSR,&firstNewSR);
		  }
	      }
	    /* Sobresale por detras el SR */
	    else
	      {
		auxSR->End = hsps[i]->EndQ;
		auxSR->Strand = firstSR->Strand;
		auxSR->Frame = firstSR->Frame;
		auxSR->Score = MAX(firstSR->Score,hsps[i]->Score);
		/*  auxSR->Locus = strcat(firstSR->Locus,hsps[i]->Subject); */
		strcpy(hsps[i]->Subject,auxSR->Locus);
		auxSR->next = NULL;
		
		allSr->nr[index]++;
		
		Chain(&newSR,auxSR,&firstNewSR);
		
		auxSR = RequestMemoryNewSR();
		
		auxSR->Start = hsps[i]->EndQ +1;
		auxSR->End = firstSR->End;
		auxSR->Strand = firstSR->Strand;
		auxSR->Frame = firstSR->Frame;
		auxSR->Score = firstSR->Score;
		strcpy(auxSR->Locus,firstSR->Locus);
		auxSR->next = NULL;

		allSr->nr[index]++;

		Chain(&newSR,auxSR,&firstNewSR);
	      }

	    /* Recompute SR chain using saved before backup pointers */
	    /*  printMess("Cut and Paste!"); */
	    numElems = CutAndPaste(&allSr->nr0[index], 
				   &firstSR,
				   &firstNewSR, 
				   &newSR,
				   &allSr->sr[index]);
	    
	    /*  printf("%d erased elems\n",numElems); */

	    allSr->nr[index] = allSr->nr[index] - numElems;
	    
	  }/* end working with sublist*/
      
      /*  printSR(allSr->sr[index],i); */
    }/* endfor */
}

void ProjectHSP(packHSP *allHsp, packSR *allSr)
{
  int i;
  
  /* projecting all the strands and frames */ 
  for(i=0; i<STRANDS*FRAMES; i++)
    project1HSP(allHsp->hsps[i],allHsp->nHsps[i],allSr);
}

#line 1458 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
void Merge(sr_t* A, sr_t* B)
{
  /* A will remain and B will be free */
  A->End = B-> End;
  A->next = B->next;

  free(B);
}

void Join1SR(packSR* allSr, int i)
{
  sr_t* A;
  sr_t* B; 
  int stop;

  A = allSr->sr[i];
  
  if (A != NULL)
    {
      stop = (A->next == NULL);
      while(!stop)
	{
	  B = A->next;
	  if (((A->End + 1) == B->Start) && (A->Score == B->Score))
	    {
	      /*  printf("Join %ld %ld -- %ld %ld\n", */
	      /*  		     A->Start,A->End, */
	      /*  		     B->Start,B->End); */
	      Merge(A,B);
	      allSr->nr[i]--;
	    }
	  else
	    {
	      /* Next loop */
	      A = A->next;
	    }
	  stop = (A->next == NULL);
	}
    }
}

void JoinSR(packSR* allSr)
{
  int i;

  /* join consecutive SR with the same score */
  for(i=0; i<STRANDS*FRAMES; i++)
    Join1SR(allSr,i);  
}

#line 1511 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
extern int VRB, GFFOUT, PSR;

extern account *m;    

/* Printing messages (information) */
void printMess(char* s)
{
  if (VRB)
     fprintf (stderr, "> %s\n",s);
}

/* Printing partial and final results (information) */
void printRes(char* s)
{
  if (VRB)
     fprintf (stderr, "\t%s\n",s);
}

/* Printing error messages */
void printError(char *s)
{
     fprintf(stderr,"Error: %s\n",s);
	  exit(1);
}

/* Get time of execution using account information */
void OutputTime()
{
  time_t tEnd;
  int t;
  float caux;
  char mess[MAXSTRING];

  (void) time(&tEnd);
  caux = (float)clock() / (float)CLOCKS_PER_SEC;

  t = (int) tEnd - m->tStart;

  if (t < caux)
    t++;

  printRes("_______________________________________________________________\n");
 
  sprintf(mess,"CPU time: \t%.3f secs",caux);
  printRes(mess);

  sprintf(mess,"Total time: \t%d secs(%.2f mins)\n\n",t,(float)t/MINUTE);
  printRes(mess);
}

void PrintHSP(packHSP* allHSP, char Query[LOCUSLENGTH])
{
  int i;
  long j;

  for(i=0; i < FRAMES*STRANDS; i++)
    for(j=0; j < allHSP->nHsps[i]; j++)
      printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\t%s\n",
	     Query,
	     VERSION,
	     HSP,
	     allHSP->hsps[i][j]->StartQ,
	     allHSP->hsps[i][j]->EndQ,
	     allHSP->hsps[i][j]->Score,
	     allHSP->hsps[i][j]->StrandQ,
	     allHSP->hsps[i][j]->FrameQ,
	     allHSP->hsps[i][j]->Subject);
}

void PrintSRStrands(packSR* allSR)
{
  int i;
  sr_t* tmp;

  printf("## blast2gff v 1.0. imim.es\n");

  for(i=0; i< FRAMES*STRANDS; i++)
    {
      printf("## Created %d SR(frame %d)\n",allSR->nr[i],i);

      tmp = allSR->sr[i];
      while(tmp != NULL)
	{
	  printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\n",
		 NONAME,
		 VERSION,
		 SR,
		 tmp->Start,
		 tmp->End,
		 tmp->Score,
		 tmp->Strand,
		 tmp->Frame); 
	  tmp = tmp->next;
	}
    }
}

void PrintSR(packSR* allSR, char Query[LOCUSLENGTH])
{
  int i;

  printf("## blast2gff v 1.0. imim.es\n");
  printf("## Created %ld SR\n",allSR->nSR);

  for(i=0; i < allSR->nSR; i++)
    {
      printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\n",
	     Query,
	     VERSION,
	     SR,
	     allSR->sortSR[i]->Start,
	     allSR->sortSR[i]->End,
	     allSR->sortSR[i]->Score,
	     allSR->sortSR[i]->Strand,
	     allSR->sortSR[i]->Frame); 
    }
}


void Output(packHSP* allHSP, packSR* allSR, char Query[LOCUSLENGTH])
{
  if (GFFOUT)
    PrintHSP(allHSP,Query);

  if (PSR)
      PrintSR(allSR,Query);
}

