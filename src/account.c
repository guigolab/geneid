/*************************************************************************
*                                                                        *
*   Module: account                                                      *
*                                                                        *
*   Accounting: results and run-time stats                               *
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

/*  $Id: account.c,v 1.5 2006-05-25 14:35:01 talioto Exp $  */

#include "geneid.h"

/* Init accounting (data structure) */
account* InitAcc()
{
  account* m; 

  /* Request memory to allocate acc. data */
  m = (account*) RequestMemoryAccounting();

  printMess("Reset accounting data");

  /* Reset counters */
  m->starts = 0;
  m->stops = 0;
  m->acc = 0;
  m->don = 0;
  m->starts_r = 0;
  m->stops_r = 0;
  m->acc_r = 0;
  m->don_r = 0;
  m->first = 0;
  m->internal = 0;
  m->terminal = 0;
  m->single = 0;
  m->orf = 0;
  m->first_r = 0;
  m->internal_r = 0;
  m->terminal_r = 0;
  m->single_r = 0;
  m->orf_r = 0;

  m->totalExons = 0;

  /* Computing the time used by every module (benchmarking): NOT USED */
  m->tSites = 0;
  m->tExons = 0;
  m->tGenes = 0;
  m->tSort = 0;
  m->tScore = 0;
  m->tBackup = 0;
  
  /* Starting the count (running time) */
  /* Real time */
  (void) time(&m->tStart);
  /* CPU time */
  clock();


  return(m);
}

/* Updating total number of sites and exons so far found in the sequence */
void updateTotals(account *m,
                  packSites* allSites,
                  packSites* allSites_r,
                  packExons* allExons,
                  packExons* allExons_r)
{
  m->starts += allSites->nStartCodons;
  m->stops += allSites->nStopCodons;
  m->acc += allSites->nAcceptorSites;
  m->acc += allSites->nU12gtagAcceptorSites;
  m->acc += allSites->nU12atacAcceptorSites;
  m->don += allSites->nDonorSites;
  m->don += allSites->nU12gtagDonorSites;
  m->don += allSites->nU12atacDonorSites;

  m->starts_r += allSites_r->nStartCodons;
  m->stops_r += allSites_r->nStopCodons;
  m->acc_r += allSites_r->nAcceptorSites;
  m->acc_r += allSites_r->nU12gtagAcceptorSites;
  m->acc_r += allSites_r->nU12atacAcceptorSites;
  m->don_r += allSites_r->nDonorSites;
  m->don += allSites_r->nU12gtagDonorSites;
  m->don += allSites_r->nU12atacDonorSites;
       
  m->first += allExons->nInitialExons;
  m->first += allExons->nU12gtagInitialExons;
  m->first += allExons->nU12atacInitialExons;
  m->internal += allExons->nInternalExons;
  m->internal += allExons->nU12gtag_U2_InternalExons;
  m->internal += allExons->nU2_U12gtag_InternalExons;
  m->internal += allExons->nU12gtag_U12gtag_InternalExons;
  m->internal += allExons->nU12atac_U2_InternalExons;
  m->internal += allExons->nU2_U12atac_InternalExons;
  m->internal += allExons->nU12atac_U12atac_InternalExons;
  m->internal += allExons->nU12gtag_U12atac_InternalExons;
  m->internal += allExons->nU12atac_U12gtag_InternalExons;
  m->terminal += allExons->nTerminalExons;
  m->terminal += allExons->nU12gtagTerminalExons;
  m->terminal += allExons->nU12atacTerminalExons;
  m->single += allExons->nSingles;
  m->orf += allExons->nORFs;

  m->first_r += allExons_r->nInitialExons;
  m->first_r += allExons_r->nU12gtagInitialExons;
  m->first_r += allExons_r->nU12atacInitialExons;
  m->internal_r += allExons_r->nInternalExons;
  m->internal_r += allExons_r->nU12gtag_U2_InternalExons;
  m->internal_r += allExons_r->nU2_U12gtag_InternalExons;
  m->internal_r += allExons_r->nU12gtag_U12gtag_InternalExons;
  m->internal_r += allExons_r->nU12atac_U2_InternalExons;
  m->internal_r += allExons_r->nU2_U12atac_InternalExons;
  m->internal_r += allExons_r->nU12atac_U12atac_InternalExons;
  m->internal_r += allExons_r->nU12gtag_U12atac_InternalExons;
  m->internal_r += allExons_r->nU12atac_U12gtag_InternalExons;
  m->terminal_r += allExons_r->nTerminalExons;
  m->terminal_r += allExons_r->nU12gtagTerminalExons;
  m->terminal_r += allExons_r->nU12atacTerminalExons;
  m->single_r += allExons_r->nSingles;  
  m->orf_r += allExons_r->nORFs;  

  m->totalExons = 
    m->first + 
    m->first_r +
    m->internal + 
    m->internal_r +
    m->terminal +
    m->terminal_r +
    m->single +
    m->single_r +
    m->orf +
    m->orf_r;
}

/* Reset acc. counters for the next input sequence */
void cleanAcc(account* m)
{
  
  /* Reset */
  m->starts = 0;
  m->stops = 0;
  m->acc = 0;
  m->don = 0;
  m->starts_r = 0;
  m->stops_r = 0;
  m->acc_r = 0;
  m->don_r = 0;

  m->first = 0;
  m->internal = 0;
  m->terminal = 0;
  m->single = 0;
  m->orf = 0;
  m->first_r = 0;
  m->internal_r = 0;
  m->terminal_r = 0;
  m->single_r = 0;
  m->orf_r = 0;

  m->totalExons = 0;
}

