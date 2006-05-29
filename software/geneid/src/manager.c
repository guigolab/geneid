/*************************************************************************
*                                                                        *
*   Module: manager                                                      *
*                                                                        *
*   Management of prediction actions: signals, exons and scores          *
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

/* $Id: manager.c,v 1.6 2006-05-29 13:47:51 talioto Exp $ */

#include "geneid.h"

extern int scanORF;
extern int U12GTAG;
extern int U12ATAC;
extern int U2GCAG;
extern long NUMSITES,NUMU12EXONS,NUMU12U12EXONS,NUMU12SITES,NUMEXONS;

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
			  packGC* GCInfo)
{
  char mess[MAXSTRING];
  int BUILD_U12_U12_INTERNALS = 0;
  
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

  allSites->nAcceptorSites =
    BuildAcceptors(Sequence,
				   sU2,
				   gp->AcceptorProfile,
				   gp->PolyPTractProfile,
				   gp->BranchPointProfile,
				   allSites->AcceptorSites,
				   l1a,l2a);
  
  sprintf(mess, "Acceptor Sites \t\t%8ld", allSites->nAcceptorSites);
  printRes(mess);

  if (U12GTAG){ 
	  allSites->nU12gtagAcceptorSites =
    	BuildU12Acceptors(Sequence,
					   sU12gtag,
					   gp->U12gtagAcceptorProfile,
					   gp->U12BranchPointProfile,
					   gp->PolyPTractProfile,
					   allSites->U12gtagAcceptorSites,
					   l1a,l2a);

	  sprintf(mess, "U12gtag Acceptor Sites \t%8ld", allSites->nU12gtagAcceptorSites);
	  printRes(mess);
  }
  if (U12ATAC){ 
	  allSites->nU12atacAcceptorSites =
    	BuildU12Acceptors(Sequence,
				 	   sU12atac,
					   gp->U12atacAcceptorProfile,
					   gp->U12BranchPointProfile,
					   gp->PolyPTractProfile,
					   allSites->U12atacAcceptorSites,
					   l1a,l2a);

	  sprintf(mess, "U12atac Acceptor Sites \t%8ld", allSites->nU12atacAcceptorSites);
	  printRes(mess);
  }  
  long numU2donsites = 0;

  allSites->nDonorSites =
    BuildDonors(Sequence,sU2, gp->DonorProfile,allSites->DonorSites,l1b,l2b,numU2donsites,NUMSITES);
  sprintf (mess,"Donor Sites \t\t%8ld", allSites->nDonorSites);
  numU2donsites = allSites->nDonorSites;
  printRes(mess);

  if (U2GCAG){
	  allSites->nDonorSites =
    	BuildDonors(Sequence,sU2gcag, gp->U2gcagDonorProfile,allSites->DonorSites,l1b,l2b,numU2donsites,NUMSITES);
	  sprintf (mess,"U2gcag Donor Sites \t%8ld", allSites->nDonorSites - numU2donsites);
	  numU2donsites = allSites->nDonorSites;
	  printRes(mess);
  }  
  if (U12GTAG){
	  allSites->nU12gtagDonorSites =
    	BuildDonors(Sequence,sU12gtag, gp->U12gtagDonorProfile,allSites->U12gtagDonorSites,l1b,l2b,0,NUMU12SITES);
	  sprintf (mess,"U12gtag Donor Sites \t%8ld", allSites->nU12gtagDonorSites);
	  printRes(mess);
  }
  if (U12ATAC){
	  allSites->nU12atacDonorSites =
    	BuildDonors(Sequence,sU12atac, gp->U12atacDonorProfile,allSites->U12atacDonorSites,l1b,l2b,0,NUMU12SITES);
	  sprintf (mess,"U12atac Donor Sites \t%8ld", allSites->nU12atacDonorSites);
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
	allSites->nU12gtagDonorSites +
	allSites->nU12gtagAcceptorSites +
	allSites->nU12atacDonorSites +
	allSites->nU12atacAcceptorSites +
	allSites->nDonorSites +	
	allSites->nStopCodons;

  sprintf(mess,"---------\t\t%8ld", allSites->nSites);
  printRes(mess);

  /* 2. Building exons with splice sites predicted before */ 
  printMess ("Computing exons ...");   
  
  /* If consecutive U12 introns are not allowed in gene model, we will not compute them. */
  if ((getkeyDict(gp->D,"U12gtag-U12gtag-Internal+")) != NOTFOUND){BUILD_U12_U12_INTERNALS = 1;}
  if ((getkeyDict(gp->D,"U12gtag-U12atac-Internal+")) != NOTFOUND){BUILD_U12_U12_INTERNALS = 1;}
  if ((getkeyDict(gp->D,"U12atac-U12gtag-Internal+")) != NOTFOUND){BUILD_U12_U12_INTERNALS = 1;}
  if ((getkeyDict(gp->D,"U12atac-U12atac-Internal+")) != NOTFOUND){BUILD_U12_U12_INTERNALS = 1;}

  allExons->nInitialExons =
    BuildInitialExons(allSites->StartCodons,allSites->nStartCodons,
					  allSites->DonorSites,allSites->nDonorSites,
					  allSites->StopCodons,allSites->nStopCodons,
					  gp->MaxDonors,sFIRST,Sequence,
					  allExons->InitialExons,NUMEXONS);
  sprintf(mess,"Initial Exons \t\t%8ld", allExons->nInitialExons);
  printRes(mess); 
  
  if (U12GTAG){  
	  allExons->nU12gtagInitialExons =
    	BuildInitialExons(allSites->StartCodons,allSites->nStartCodons,
						  allSites->U12gtagDonorSites,allSites->nU12gtagDonorSites,
						  allSites->StopCodons,allSites->nStopCodons,
						  gp->MaxDonors,sU12gtagFIRST,Sequence,
						  allExons->U12gtagInitialExons,NUMU12EXONS);
	  sprintf(mess,"U12gtag Initial Exons \t%8ld", allExons->nU12gtagInitialExons);
	  printRes(mess); 
  }
   if (U12ATAC){  
	  allExons->nU12atacInitialExons =
    	BuildInitialExons(allSites->StartCodons,allSites->nStartCodons,
						  allSites->U12atacDonorSites,allSites->nU12atacDonorSites,
						  allSites->StopCodons,allSites->nStopCodons,
						  gp->MaxDonors,sU12atacFIRST,Sequence,
						  allExons->U12atacInitialExons,NUMU12EXONS);
	  sprintf(mess,"U12atac Initial Exons \t%8ld", allExons->nU12atacInitialExons);
	  printRes(mess); 
  } 
  allExons->nInternalExons =
	BuildInternalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
					   allSites->DonorSites,allSites->nDonorSites,
					   allSites->StopCodons,allSites->nStopCodons,
					   gp->MaxDonors,sINTERNAL,Sequence,
					   allExons->InternalExons,NUMEXONS);
  sprintf(mess,"Internal Exons \t\t%8ld", allExons->nInternalExons);
  printRes(mess); 

  if (U12GTAG){   
	  if (BUILD_U12_U12_INTERNALS){
	  	allExons->nU12gtag_U12gtag_InternalExons =
		BuildInternalExons(allSites->U12gtagAcceptorSites,allSites->nU12gtagAcceptorSites,
						   allSites->U12gtagDonorSites,allSites->nU12gtagDonorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   gp->MaxDonors,sU12gtag_U12gtagINTERNAL,Sequence,
						   allExons->U12gtag_U12gtag_InternalExons,NUMU12EXONS);
	  	sprintf(mess,"U12gtag-U12gtag Internal Exons \t%8ld", allExons->nU12gtag_U12gtag_InternalExons);
	  	printRes(mess); 
	  }
	  allExons->nU12gtag_U2_InternalExons =
		BuildInternalExons(allSites->U12gtagAcceptorSites,allSites->nU12gtagAcceptorSites,
						   allSites->DonorSites,allSites->nDonorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   gp->MaxDonors,sU12gtag_U2INTERNAL,Sequence,
						   allExons->U12gtag_U2_InternalExons,NUMU12EXONS);
	  sprintf(mess,"U12gtag-U2 Internal Exons \t%8ld", allExons->nU12gtag_U2_InternalExons);
	  printRes(mess); 

	  allExons->nU2_U12gtag_InternalExons =
		BuildInternalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
						   allSites->U12gtagDonorSites,allSites->nU12gtagDonorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   gp->MaxDonors,sU2_U12gtagINTERNAL,Sequence,
						   allExons->U2_U12gtag_InternalExons,NUMU12EXONS);
	  sprintf(mess,"U2-U12gtag Internal Exons \t%8ld", allExons->nU2_U12gtag_InternalExons);
	  printRes(mess); 
  }
  if (U12ATAC){   
   	  if (BUILD_U12_U12_INTERNALS){
	  	allExons->nU12atac_U12atac_InternalExons =
		BuildInternalExons(allSites->U12atacAcceptorSites,allSites->nU12atacAcceptorSites,
						   allSites->U12atacDonorSites,allSites->nU12atacDonorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   gp->MaxDonors,sU12atac_U12atacINTERNAL,Sequence,
						   allExons->U12atac_U12atac_InternalExons,NUMU12EXONS);
	  	sprintf(mess,"U12atac-U12atac Internal Exons \t%8ld", allExons->nU12atac_U12atac_InternalExons);
	  	printRes(mess); 
	  }
	  allExons->nU12atac_U2_InternalExons =
		BuildInternalExons(allSites->U12atacAcceptorSites,allSites->nU12atacAcceptorSites,
						   allSites->DonorSites,allSites->nDonorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   gp->MaxDonors,sU12atac_U2INTERNAL,Sequence,
						   allExons->U12atac_U2_InternalExons,NUMU12EXONS);
	  sprintf(mess,"U12atac-U2 Internal Exons \t%8ld", allExons->nU12atac_U2_InternalExons);
	  printRes(mess); 

	  allExons->nU2_U12atac_InternalExons =
		BuildInternalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
						   allSites->U12atacDonorSites,allSites->nU12atacDonorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   gp->MaxDonors,sU2_U12atacINTERNAL,Sequence,
						   allExons->U2_U12atac_InternalExons,NUMU12EXONS);
	  sprintf(mess,"U2-U12atac Internal Exons \t%8ld", allExons->nU2_U12atac_InternalExons);
	  printRes(mess); 
  }
  if (U12GTAG && U12ATAC){  
	  if (BUILD_U12_U12_INTERNALS){ 
	  	allExons->nU12gtag_U12atac_InternalExons =
		BuildInternalExons(allSites->U12gtagAcceptorSites,allSites->nU12gtagAcceptorSites,
						   allSites->U12atacDonorSites,allSites->nU12atacDonorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   gp->MaxDonors,sU12gtag_U12atacINTERNAL,Sequence,
						   allExons->U12gtag_U12atac_InternalExons,NUMU12U12EXONS);
	  	sprintf(mess,"U12gtag-U12atac Internal Exons \t%8ld", allExons->nU12gtag_U12atac_InternalExons);
	  	printRes(mess); 

	  	allExons->nU12atac_U12gtag_InternalExons =
		BuildInternalExons(allSites->U12atacAcceptorSites,allSites->nU12atacAcceptorSites,
						   allSites->U12gtagDonorSites,allSites->nU12gtagDonorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   gp->MaxDonors,sU12atac_U12gtagINTERNAL,Sequence,
						   allExons->U12atac_U12gtag_InternalExons,NUMU12U12EXONS);
	  	sprintf(mess,"U12atac-U12gtag Internal Exons \t%8ld", allExons->nU12atac_U12gtag_InternalExons);
	  	printRes(mess); 
	  }
  } 
  
  allExons->nTerminalExons =
    BuildTerminalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
					   allSites->StopCodons,allSites->nStopCodons,
					   LengthSequence,cutPoint,sTERMINAL,Sequence,
					   allExons->TerminalExons,NUMEXONS);
  sprintf(mess,"Terminal Exons \t\t%8ld", allExons->nTerminalExons);
  printRes(mess); 

  if (U12GTAG){ 
	  allExons->nU12gtagTerminalExons =
    	BuildTerminalExons(allSites->U12gtagAcceptorSites,allSites->nU12gtagAcceptorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   LengthSequence,cutPoint,sU12gtagTERMINAL,Sequence,
						   allExons->U12gtagTerminalExons,NUMU12EXONS);
	  sprintf(mess,"U12gtag Terminal Exons \t%8ld", allExons->nU12gtagTerminalExons);
	  printRes(mess); 
  }

  if (U12ATAC){ 
	  allExons->nU12atacTerminalExons =
    	BuildTerminalExons(allSites->U12atacAcceptorSites,allSites->nU12atacAcceptorSites,
						   allSites->StopCodons,allSites->nStopCodons,
						   LengthSequence,cutPoint,sU12atacTERMINAL,Sequence,
						   allExons->U12atacTerminalExons,NUMU12EXONS);
	  sprintf(mess,"U12atac Terminal Exons \t%8ld", allExons->nU12atacTerminalExons);
	  printRes(mess); 
  }
  
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
	allExons->nU12gtagInitialExons +
	allExons->nU12atacInitialExons +    
	allExons->nU12gtag_U12gtag_InternalExons +
	allExons->nU12gtag_U2_InternalExons +
	allExons->nU2_U12gtag_InternalExons +
	allExons->nU12atac_U12atac_InternalExons +
	allExons->nU12atac_U2_InternalExons +
	allExons->nU2_U12atac_InternalExons +
    allExons->nU12gtag_U12atac_InternalExons +
	allExons->nU12atac_U12gtag_InternalExons +
	allExons->nU12gtagTerminalExons +
    allExons->nU12atacTerminalExons +
    allExons->nSingles +
    allExons->nORFs;

  sprintf(mess,"---------\t\t%8ld", allExons->nExons);
  printRes(mess); 
}
