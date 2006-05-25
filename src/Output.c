/*************************************************************************
*                                                                        *
*   Module: Output                                                       *
*                                                                        *
*   Management: displaying results                                       *
*                                                                        *
*   This file is part of the geneid 1.2 distribution                     *
*                                                                        *
*     Copyright (C) 2003 - Enrique BLANCO GARCIA                         *
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

/*  $Id: Output.c,v 1.13 2006-05-25 14:15:44 talioto Exp $  */

#include "geneid.h"
extern int U12GTAG;
extern int U12ATAC;
extern int VRB;

extern int SFP, SDP, SAP, STP,
           EFP, EIP, ETP, EXP, ESP,EOP,
           FWD, RVS, EVD,
           GENAMIC, GENEID,
           GFF, GFF3, XML;

extern account *m;

/* Printing messages (information) */
void printMess(char* s)
{
  if (VRB)
    fprintf (stderr, "> %s\n",s);
}

/* Printing stats and results (information) */
void printRes(char* s)
{
  if (VRB)
    fprintf (stderr, "\t%s\n",s);
}

/* Printing error messages */
void printError(char* s)
{
  fprintf(stderr,"Error: %s\n",s);
  exit(1);
}

/* Printing messages (nucleotides read) */
void printReadingInfo(char* s)
{
  if (VRB)
    fprintf (stderr, "%s\n",s);
}

/* Printing messages (parameters read for PWA - signal prediction) */
void PrintProfile (profile* p, char* signal)
{
  char mess[MAXSTRING];

  sprintf(mess, 
		  "Reading... %s:\t%d\t%d\t%d\t(%ld)\t%5.2f", 
		  signal,
		  p->dimension, 
		  p->offset, 
		  p->order,
		  p->dimensionTrans,
		  p->cutoff);
  
  printMess(mess);
}

/* Output: header for results displayed immediately */
void OutputHeader(char* locus, long l)
{
  char* s;
  char mess[MAXSTRING];

  /* 0. Size checkpoint and information */
  if (!l)
     {
	   sprintf(mess,"%s: sequence is empty",locus);     
	   printError(mess);
     }
  else
	{
	  sprintf(mess,"%s: %ld nucleotides\n",locus,l);     
	  printMess(mess);
	}
  
  /* 1. Extract the starting time to display */
  s = ctime(&m->tStart);
  
  /* 2. Output headers: gff, geneid or xml format */
  
  if (GFF3){
    printf("## gff-version 3\n");
  } else {
  	if (GFF)
    printf("## gff-version 2\n");
  }
  
  if (XML)
    {
      /* XML format header */
      printf("<?xml version=\"1.0\" ?>\n");  
      printf("<!DOCTYPE prediction SYSTEM \"geneid.dtd\">\n");
      s[strlen(s)-1] = '\0';  
      printf("<prediction locus=\"%s\" length=\"%ld\" source=\"%s\" date=\"%s\"",
			 locus,l,VERSION,s);
    }   
  else
    {
      /* gff and geneid formats */
      s[strlen(s)-1] = '\n';
      if (GFF3){  
	  	printf("# date %s",s);
      	printf("# source-version: %s -- geneid@imim.es\n",VERSION);
      	printf("##sequence-region %s 1 %ld\n",locus,l);
	  } else {
	  	printf("## date %s",s);
      	printf("## source-version: %s -- geneid@imim.es\n",VERSION);
      	printf("# Sequence %s - Length = %ld bps\n",locus,l);
	  }
    }
}

/* Display some predictions results according to the options selected */
void Output(packSites* allSites,
            packSites* allSites_r,
            packExons* allExons,
            packExons* allExons_r,
            exonGFF* exons,
            long nExons,
            char* Locus,
            long l1,
            long l2,
            char* Sequence,
            gparam* gp,
            dict* dAA)
{
  /* 1. Printing Forward */
  if (FWD)
    {
      printMess("Printing forward selected elements");
	  
      /* sites */
      if (SFP) 
		PrintSites(allSites->StartCodons, allSites->nStartCodons,
				   STA, Locus, FORWARD, l1, l2, Sequence, gp->StartProfile);
      if (SAP){
		PrintSites(allSites->AcceptorSites, allSites->nAcceptorSites,
				   ACC, Locus, FORWARD, l1, l2, Sequence, gp->AcceptorProfile);
		if (U12GTAG){
			PrintSites(allSites->U12gtagAcceptorSites, allSites->nU12gtagAcceptorSites,
				   ACC, Locus, FORWARD, l1, l2, Sequence, gp->U12gtagAcceptorProfile);
		}		   
		if (U12ATAC){
			PrintSites(allSites->U12atacAcceptorSites, allSites->nU12atacAcceptorSites,
				   ACC, Locus, FORWARD, l1, l2, Sequence, gp->U12atacAcceptorProfile);
		}
	  }
      if (SDP){
		PrintSites(allSites->DonorSites, allSites->nDonorSites,
				   DON, Locus, FORWARD, l1, l2, Sequence, gp->DonorProfile);
		if (U12GTAG){
			PrintSites(allSites->U12gtagDonorSites, allSites->nU12gtagDonorSites,
				   DON, Locus, FORWARD, l1, l2, Sequence, gp->U12gtagDonorProfile);
		}	
		if (U12ATAC){
			PrintSites(allSites->U12atacDonorSites, allSites->nU12atacDonorSites,
				   DON, Locus, FORWARD, l1, l2, Sequence, gp->U12atacDonorProfile);
		}	   
	  }
      if (STP)
		PrintSites(allSites->StopCodons, allSites->nStopCodons,
				   STO, Locus, FORWARD, l1, l2, Sequence, gp->StopProfile);
      
      /* exons */
      if (EFP){
		PrintExons(allExons->InitialExons,allExons->nInitialExons,
				   FIRST, Locus, l1, l2, Sequence, dAA);
		if (U12GTAG){
			PrintExons(allExons->U12gtagInitialExons,allExons->nU12gtagInitialExons,
				   FIRST, Locus, l1, l2, Sequence, dAA);
		}
		if (U12ATAC){
			PrintExons(allExons->U12atacInitialExons,allExons->nU12atacInitialExons,
				   FIRST, Locus, l1, l2, Sequence, dAA);
		}
      }
	  if (EIP){
		PrintExons(allExons->InternalExons,allExons->nInternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
		if (U12GTAG){
			PrintExons(allExons->U2_U12gtag_InternalExons,allExons->nU2_U12gtag_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons->U12gtag_U2_InternalExons,allExons->nU12gtag_U2_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons->U12gtag_U12gtag_InternalExons,allExons->nU12gtag_U12gtag_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
		}
		if (U12ATAC){
			PrintExons(allExons->U2_U12atac_InternalExons,allExons->nU2_U12atac_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons->U12atac_U2_InternalExons,allExons->nU12atac_U2_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons->U12atac_U12atac_InternalExons,allExons->nU12atac_U12atac_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
		}
		if (U12ATAC && U12GTAG){
			PrintExons(allExons->U12gtag_U12atac_InternalExons,allExons->nU12gtag_U12atac_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons->U12atac_U12gtag_InternalExons,allExons->nU12atac_U12gtag_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
		}
      }
	  if (ETP){
		PrintExons(allExons->TerminalExons,allExons->nTerminalExons,
				   TERMINAL, Locus, l1, l2, Sequence, dAA);
		if (U12GTAG){
				PrintExons(allExons->U12gtagTerminalExons,allExons->nU12gtagTerminalExons,
				   TERMINAL, Locus, l1, l2, Sequence, dAA);
		}
		if (U12ATAC){
				PrintExons(allExons->U12atacTerminalExons,allExons->nU12atacTerminalExons,
				   TERMINAL, Locus, l1, l2, Sequence, dAA);
		}
      }
	  if (ESP)
		PrintExons(allExons->Singles,allExons->nSingles,
				   SINGLE, Locus, l1, l2, Sequence, dAA);
      if (EOP)
		PrintExons(allExons->ORFs,allExons->nORFs,
				   ORF, Locus, l1, l2, Sequence, dAA);
    }
  
  /* 2. Printing Reverse */
  if (RVS)
    {
      printMess("Printing reverse selected elements\n");
	  
      /* sites */ 
      if (SFP)
		PrintSites(allSites_r->StartCodons,allSites_r->nStartCodons,STA,
				   Locus,REVERSE, l1, l2, Sequence, gp->StartProfile);
      if (SAP){
		PrintSites(allSites_r->AcceptorSites, allSites_r->nAcceptorSites,
				   ACC, Locus, REVERSE, l1, l2, Sequence, gp->AcceptorProfile);
		if (U12GTAG){
			PrintSites(allSites_r->U12gtagAcceptorSites, allSites_r->nU12gtagAcceptorSites,
				   ACC, Locus, REVERSE, l1, l2, Sequence, gp->U12gtagAcceptorProfile);
		}		   
		if (U12ATAC){
			PrintSites(allSites_r->U12atacAcceptorSites, allSites_r->nU12atacAcceptorSites,
				   ACC, Locus, REVERSE, l1, l2, Sequence, gp->U12atacAcceptorProfile);
		}
	  }
      if (SDP){
		PrintSites(allSites_r->DonorSites, allSites_r->nDonorSites,
				   DON, Locus, REVERSE, l1, l2, Sequence, gp->DonorProfile);
		if (U12GTAG){
			PrintSites(allSites_r->U12gtagDonorSites, allSites_r->nU12gtagDonorSites,
				   DON, Locus, REVERSE, l1, l2, Sequence, gp->U12gtagDonorProfile);
		}	
		if (U12ATAC){
			PrintSites(allSites_r->U12atacDonorSites, allSites_r->nU12atacDonorSites,
				   DON, Locus, REVERSE, l1, l2, Sequence, gp->U12atacDonorProfile);
		}	   
	  }
      if (STP)
		PrintSites(allSites_r->StopCodons,allSites_r->nStopCodons,STO,
				   Locus,REVERSE, l1, l2, Sequence, gp->StopProfile);
      
      /* exons */
      if (EFP){
		PrintExons(allExons_r->InitialExons,allExons_r->nInitialExons,
				   FIRST, Locus, l1, l2, Sequence, dAA);
		if (U12GTAG){
			PrintExons(allExons_r->U12gtagInitialExons,allExons_r->nU12gtagInitialExons,
				   FIRST, Locus, l1, l2, Sequence, dAA);
		}
		if (U12ATAC){
			PrintExons(allExons_r->U12atacInitialExons,allExons_r->nU12atacInitialExons,
				   FIRST, Locus, l1, l2, Sequence, dAA);
		}
      }
	  if (EIP){
		PrintExons(allExons_r->InternalExons,allExons_r->nInternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
		if (U12GTAG){
			PrintExons(allExons_r->U2_U12gtag_InternalExons,allExons_r->nU2_U12gtag_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons_r->U12gtag_U2_InternalExons,allExons_r->nU12gtag_U2_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons_r->U12gtag_U12gtag_InternalExons,allExons_r->nU12gtag_U12gtag_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
		}
		if (U12ATAC){
			PrintExons(allExons_r->U2_U12atac_InternalExons,allExons_r->nU2_U12atac_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons_r->U12atac_U2_InternalExons,allExons_r->nU12atac_U2_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons_r->U12atac_U12atac_InternalExons,allExons_r->nU12atac_U12atac_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
		}
		if (U12ATAC && U12GTAG){
			PrintExons(allExons_r->U12gtag_U12atac_InternalExons,allExons_r->nU12gtag_U12atac_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
			PrintExons(allExons_r->U12atac_U12gtag_InternalExons,allExons_r->nU12atac_U12gtag_InternalExons,
				   INTERNAL, Locus, l1, l2, Sequence, dAA);
		}
      }
	  if (ETP){
		PrintExons(allExons_r->TerminalExons,allExons_r->nTerminalExons,
				   TERMINAL, Locus, l1, l2, Sequence, dAA);
		if (U12GTAG){
				PrintExons(allExons_r->U12gtagTerminalExons,allExons_r->nU12gtagTerminalExons,
				   TERMINAL, Locus, l1, l2, Sequence, dAA);
		}
		if (U12ATAC){
				PrintExons(allExons_r->U12atacTerminalExons,allExons_r->nU12atacTerminalExons,
				   TERMINAL, Locus, l1, l2, Sequence, dAA);
		}
      }
      if (ESP)
		PrintExons(allExons_r->Singles,allExons_r->nSingles,
				   SINGLE, Locus, l1, l2, Sequence, dAA);
      if (EOP)
		PrintExons(allExons_r->ORFs,allExons_r->nORFs,
				   ORF, Locus, l1, l2, Sequence, dAA);
    }
  
  /* 3. Print all exons */
  if (EXP) 
    {
      printMess("Printing all predicted Exons of current split\n");   
      PrintExons(exons, nExons, FIRST + INTERNAL + TERMINAL + SINGLE + ORF, 
				 Locus, l1, l2, Sequence, dAA);
    }
}

/* Print best genes using selected format */
void OutputGene(packGenes* pg,
                long nExons,
                char* Locus,
                char* Sequence,
                gparam* gp,
                dict* dAA)
{
  /* Retrieving the best predicted genes recursively */
  if (nExons>0)
    {
      printMess("Recovering gene-solution...");
      CookingGenes(pg->GOptim, Locus, Sequence, gp, dAA);
      if (XML)
		printf("</prediction>\n");
    }
  else
	if (XML)
	  printf(" genes=\"0\" score =\"0.00\">\n</prediction>\n");	
}

/* Display information about stats of predicted sites and exons */
void OutputStats(char* Locus)
{
  char mess[MAXSTRING];
  
  if (GENEID)
    {
      sprintf(mess,"\n\tStats (Sequence %s)",Locus);
      printRes(mess);
      
      printRes("__________________________________________________________________________\n");
	  
      sprintf(mess,"%8s\t%8s\t%8s\t%8s\t%8s",
			  sFIRST,sINTERNAL,sTERMINAL,sSINGLE,sORF);
      
      printRes(mess);
      printRes("__________________________________________________________________________\n");
	  
      sprintf(mess,"%8ld\t%8ld\t%8ld\t%8ld\t%8ld\n",      
			  m->first, 
			  m->internal, 
			  m->terminal,
			  m->single,
			  m->orf); 
      printRes(mess);
	  
      sprintf(mess,"%8ld\t%8ld\t%8ld\t%8ld\t%8ld\n",      
			  m->first_r, 
			  m->internal_r, 
			  m->terminal_r,
			  m->single_r,
			  m->orf_r);      
      printRes(mess);
	  
      sprintf(mess,"TOTAL: %ld predicted exons\n",
			  m->totalExons);
      
      printRes(mess);
    }
}

/* Computing running time by using accounting information */
void OutputTime()
{
  time_t tEnd;
  int t;
  float caux;
  char mess[MAXSTRING];

  /* Final time */
  /* Real time */
  (void) time(&tEnd);
  /* CPU time */
  caux = (float)clock() / (float)CLOCKS_PER_SEC;

  t = (int) tEnd - m->tStart;

  /* Correction */
  if (t < caux)
    t++;

  printRes("__________________________________________________________________________\n");
  
  sprintf(mess,"CPU time: \t%.3f secs",caux);
  printRes(mess);
  
  sprintf(mess,"Total time: \t%d secs(%.2f mins)\n\n",t,(float)t/MINUTE);
  printRes(mess);
}



