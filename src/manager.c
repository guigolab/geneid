/*************************************************************************
*                                                                        *
*   Module: manager                                                      *
*                                                                        *
*   Management of prediction actions: signals, exons and scores          *
*                                                                        *
*   This file is part of the geneid 1.3 distribution                     *
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

/* $Id: manager.c,v 1.10 2006-12-21 13:56:54 talioto Exp $ */

#include "geneid.h"

extern int scanORF;
extern int U12GTAG;
extern int U12ATAC;
extern int U2GCAG;
extern int U2GTA;
extern int U2GTG;
extern int U2GTY;
extern long NUMSITES,NUMEXONS;

/* Management of splice sites prediction and exon construction/scoring */
void  manager(char *Sequence, 
	      long LengthSequence,
	      packSites* allSites,
	      packExons* allExons,
	      long l1, long l2, long lowerlimit, long upperlimit,
	      int Strand,
	      packExternalInformation* external,
	      packHSP* hsp,
	      gparam* gp,
	      gparam** isochores,
	      int nIsochores,
	      packGC* GCInfo,
	      site* acceptorsites,
	      site* donorsites
	      )
{
  char mess[MAXSTRING];

  /* For sorting sites */
/*   site* acceptorsites;  */
/*   site* donorsites;  */
  long l1a, l1b,
	l2a, l2b,
	l1c, l2c;

  long cutPoint;
 
  /* 0. Define boundaries of splice site prediction
	 according to current split positions and strand selected */
  if (Strand == FORWARD)
    {
      /* Forward sense */
      /* Start codons and Acceptor sites limits */
      l1a = l1;
      l2a = (l2 == upperlimit)? l2 : l2 - OVERLAP;

      /* Donor sites limits */
      l1b = l1;
      l2b = l2;

      /* Stop codon limits */
      l1c = l1;
      l2c = l2;

      /* Terminal/Single exons: */
      /* are allowed if their Stop codon is placed behind cutPoint */
      /* FWD: every stop codon might be used without problems */
      cutPoint = l1;
    }
  else
    {
      /* Reverse sense */
      /* Start codons and Acceptor sites limits */
      l1a = l1;
      l2a = l2;

      /* Donor sites limits */
      l1b = (l1 == lowerlimit)? l1: l1 + OVERLAP;
      l2b = l2;

      /* Stop codon limits */
      l1c = l1;
      l2c = l2;

      /* Terminal/Single exons: */
      /* are allowed if their Stop codon is placed behind cutPoint (RVS) */
      /* RVS: reading from right to left the forward sense sequence */
      cutPoint = (l1 == lowerlimit)? l1 : l1 + OVERLAP;
    }

  /* 1. Predicting splice sites of current split of DNA sequence */ 
  printMess ("Computing sites ...");

  allSites->nStartCodons =
    GetSitesWithProfile(Sequence,gp->StartProfile,allSites->StartCodons,l1a,l2a);
  sprintf(mess, "Start Codons \t\t%8ld", allSites->nStartCodons);
  printRes(mess);
  
  long numAccsites = 0;
  
  allSites->nAcceptorSites =
    BuildAcceptors(Sequence,
		   U2,
		   sU2type,
		   sU2,
		   gp->AcceptorProfile,
		   gp->PolyPTractProfile,
		   gp->BranchPointProfile,
		   allSites->AcceptorSites,
		   l1a,l2a,numAccsites,NUMSITES);
  
  sprintf(mess, "Acceptor Sites \t\t%8ld", allSites->nAcceptorSites - numAccsites);
  numAccsites = allSites->nAcceptorSites;
  printRes(mess);

  if (U12GTAG){ 
	  allSites->nAcceptorSites =
	    BuildU12Acceptors(Sequence,U12gtag,sU12type,
					   sU12gtag,
					   gp->U12gtagAcceptorProfile,
					   gp->U12BranchPointProfile,
					   gp->PolyPTractProfile,
					   allSites->AcceptorSites,
					   l1a,l2a,numAccsites,NUMSITES);

	  sprintf(mess, "U12gtag Acceptor Sites \t%8ld", allSites->nAcceptorSites - numAccsites);
	  numAccsites = allSites->nAcceptorSites;
	  printRes(mess);
  }
  if (U12ATAC){ 
	  allSites->nAcceptorSites =
	    BuildU12Acceptors(Sequence,U12atac,sU12type,
				 	   sU12atac,
					   gp->U12atacAcceptorProfile,
					   gp->U12BranchPointProfile,
					   gp->PolyPTractProfile,
					   allSites->AcceptorSites,
					   l1a,l2a,numAccsites,NUMSITES);

	  sprintf(mess, "U12atac Acceptor Sites \t%8ld", allSites->nAcceptorSites - numAccsites);
	  numAccsites = allSites->nAcceptorSites;
	  printRes(mess);
  }  

  long numDonsites = 0;

  allSites->nDonorSites =
    BuildDonors(Sequence,U2,sU2type,sU2, gp->DonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES);
  sprintf (mess,"Donor Sites \t\t%8ld", allSites->nDonorSites);
  numDonsites = allSites->nDonorSites;
  printRes(mess);

  if (U12GTAG){
	  allSites->nDonorSites =
	    BuildDonors(Sequence,U12gtag,sU12type,sU12gtag, gp->U12gtagDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES);
	  sprintf (mess,"U12gtag Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }
  if (U12ATAC){
	  allSites->nDonorSites =
	    BuildDonors(Sequence, U12atac,sU12type,sU12atac, gp->U12atacDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES);
	  sprintf (mess,"U12atac Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }
  if (U2GCAG){
	  allSites->nDonorSites =
	    BuildDonors(Sequence,U2, sU2type,sU2gcag, gp->U2gcagDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES);
	  sprintf (mess,"U2gcag Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }  
  if (U2GTA){
	  allSites->nDonorSites =
    	BuildDonors(Sequence,U2,sU2type,sU2gta, gp->U2gtaDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES);
	  sprintf (mess,"U2gta Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }  
  if (U2GTG){
	  allSites->nDonorSites =
    	BuildDonors(Sequence,U2,sU2type,sU2gtg, gp->U2gtgDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES);
	  sprintf (mess,"U2gtg Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }  
  if (U2GTY){
	  allSites->nDonorSites =
    	BuildDonors(Sequence,U2,sU2type,sU2gty, gp->U2gtyDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES);
	  sprintf (mess,"U2gty Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }  

  allSites->nStopCodons =
	GetStopCodons(Sequence,gp->StopProfile, allSites->StopCodons,l1c,l2c);
  sprintf (mess,"Stop Codons \t\t%8ld", allSites->nStopCodons);
  printRes(mess);
  
  /* Total number of predicted splice sites in this strand */
  allSites->nSites =
	allSites->nStartCodons +
	allSites->nAcceptorSites +
	allSites->nDonorSites +	
	allSites->nStopCodons;

  sprintf(mess,"---------\t\t%8ld", allSites->nSites);
  printRes(mess);

  /* Predicted sites must be sorted by position */
  printMess ("Sorting sites ...");
  SortSites(allSites->DonorSites,allSites->nDonorSites,donorsites,l1b,l2b);
  SortSites(allSites->AcceptorSites,allSites->nAcceptorSites,acceptorsites,l1a,l2a);

  /* 2. Building exons with splice sites predicted before */ 
  printMess ("Computing exons ...");   
  

  allExons->nInitialExons =
    BuildInitialExons(allSites->StartCodons,allSites->nStartCodons,
					  allSites->DonorSites,allSites->nDonorSites,
					  allSites->StopCodons,allSites->nStopCodons,
					  gp->MaxDonors,sFIRST,Sequence,
					  allExons->InitialExons,NUMEXONS);
  sprintf(mess,"Initial Exons \t\t%8ld", allExons->nInitialExons);
  printRes(mess); 

  allExons->nInternalExons =
	BuildInternalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
					   allSites->DonorSites,allSites->nDonorSites,
					   allSites->StopCodons,allSites->nStopCodons,
					   gp->MaxDonors,sINTERNAL,Sequence,
					   allExons->InternalExons,NUMEXONS);
  sprintf(mess,"Internal Exons \t\t%8ld", allExons->nInternalExons);
  printRes(mess); 
  
  allExons->nTerminalExons =
    BuildTerminalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
					   allSites->StopCodons,allSites->nStopCodons,
					   LengthSequence,cutPoint,sTERMINAL,Sequence,
					   allExons->TerminalExons,NUMEXONS);
  sprintf(mess,"Terminal Exons \t\t%8ld", allExons->nTerminalExons);
  printRes(mess); 
  
  allExons->nSingles =
    BuildSingles(allSites->StartCodons,allSites->nStartCodons,
				 allSites->StopCodons,allSites->nStopCodons,
				 cutPoint, Sequence,
				 allExons->Singles);
  sprintf(mess,"Single genes \t\t%8ld", allExons->nSingles);
  printRes(mess); 

  if (scanORF)
    {
      allExons->nORFs =
        BuildORFs(allSites->StopCodons,allSites->nStopCodons,
				  allSites->StopCodons,allSites->nStopCodons,
				  cutPoint, Sequence,
				  allExons->ORFs);
      sprintf(mess,"ORFs \t\t\t%8ld", allExons->nORFs);
      printRes(mess); 
    }
  else
	allExons->nORFs = 0;

  /* 3. Scoring and Filtering Exons */
  ScoreExons(Sequence, allExons, 
             l1, l2, Strand, 
	     external, hsp,
             isochores,nIsochores,
             GCInfo);
  
  /* Total number of built exons in this strand */
  allExons->nExons =
    allExons->nInitialExons +
    allExons->nInternalExons +
    allExons->nTerminalExons +
    allExons->nSingles +
    allExons->nORFs;

  sprintf(mess,"---------\t\t%8ld", allExons->nExons);
  printRes(mess); 
}
