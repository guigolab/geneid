/*************************************************************************
*                                                                        *
*   Module: manager                                                      *
*                                                                        *
*   Management of prediction actions: signals, exons and scores          *
*                                                                        *
*   This file is part of the geneid 1.1 distribution                     *
*                                                                        *
*     Copyright (C) 2001 - Enrique BLANCO GARCIA                         *
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

/* $Id: manager.c,v 1.1 2001-12-18 16:22:36 eblanco Exp $ */

#include "geneid.h"

extern int scanORF;

/* Management of splice sites prediction and exon construction/scoring */
void  manager(char *Sequence, 
			  long LengthSequence,
			  packSites* allSites,
			  packExons* allExons,
			  long l1, long l2,
			  int Strand,
			  packSR* sr,
			  gparam* gp,
			  gparam** isochores,
			  int nIsochores,
			  packGC* GCInfo)
{
  char mess[MAXSTRING];
  
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
      l2a = (l2 == LengthSequence - 1)? l2 : l2 - OVERLAP;

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
      l1b = (l1 == 0)? l1: l1 + OVERLAP;
      l2b = l2;

      /* Stop codon limits */
      l1c = l1;
      l2c = l2;

      /* Terminal/Single exons: */
      /* are allowed if their Stop codon is placed behind cutPoint (RVS) */
      /* RVS: reading from right to left the forward sense sequence */
      cutPoint = (l1 == 0)? l1 : l1 + OVERLAP;
    }

  /* 1. Predicting splice sites of current split of DNA sequence */ 
  printMess ("Computing sites ...");

  allSites->nStartCodons =
    GetSitesWithProfile(Sequence,gp->StartProfile,allSites->StartCodons,l1a,l2a);
  sprintf(mess, "Start Codons \t\t%8ld", allSites->nStartCodons);
  printRes(mess);

  allSites->nAcceptorSites =
    GetSitesWithProfile(Sequence,gp->AcceptorProfile,allSites->AcceptorSites,l1a,l2a);
  sprintf(mess, "Acceptor Sites \t\t%8ld", allSites->nAcceptorSites);
  printRes(mess);

  allSites->nDonorSites =
    GetSitesWithProfile(Sequence,gp->DonorProfile,allSites->DonorSites,l1b,l2b);
  sprintf (mess,"Donor Sites \t\t%8ld", allSites->nDonorSites);
  printRes(mess);

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

  /* 2. Building exons with splice sites predicted before */ 
  printMess ("Computing exons ...");   
  
  allExons->nInitialExons =
    BuildInitialExons(allSites->StartCodons,allSites->nStartCodons,
					  allSites->DonorSites,allSites->nDonorSites,
					  allSites->StopCodons,allSites->nStopCodons,
					  gp->MaxDonors,Sequence,
					  allExons->InitialExons);
  sprintf(mess,"Initial Exons \t\t%8ld", allExons->nInitialExons);
  printRes(mess); 
  
  allExons->nInternalExons =
	BuildInternalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
					   allSites->DonorSites,allSites->nDonorSites,
					   allSites->StopCodons,allSites->nStopCodons,
					   gp->MaxDonors,Sequence,
					   allExons->InternalExons);
  sprintf(mess,"Internal Exons \t\t%8ld", allExons->nInternalExons);
  printRes(mess); 
  
  allExons->nTerminalExons =
    BuildTerminalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
					   allSites->StopCodons,allSites->nStopCodons,
					   LengthSequence,cutPoint,Sequence,
					   allExons->TerminalExons);
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

  /* 3. Scoring and Filtering Exons */
  ScoreExons(Sequence, allExons, 
             l1, l2, Strand, sr,
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
