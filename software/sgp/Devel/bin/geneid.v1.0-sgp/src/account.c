/*************************************************************************
*                                                                        *
*   Module: account                                                      *
*                                                                        *
*   Accounting of results stats and time of execution.                   *
*                                                                        *
*   This file is part of the geneid Distribution                         *
*                                                                        *
*     Copyright (C) 2000 - Enrique BLANCO GARCIA                         *
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

/*  $Id: account.c,v 1.1 2000-07-05 08:33:26 eblanco Exp $  */

#include "geneid.h"

/* Init accounting */
account* InitAcc()
{
  account* m; 
  
  m = (account *) malloc(sizeof(account));

  /* Reset counters */
  m->totalExons = 0;
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
  m->first_r = 0;
  m->internal_r = 0;
  m->terminal_r = 0;
  m->single_r = 0;

  m->tSites = 0;
  m->tExons = 0;
  m->tGenes = 0;
  m->tSort = 0;
  m->tScore = 0;
  m->tBackup = 0;
  
  printMess("Reset accounting");
  /* What time is it? */ 
  clock();
  (void) time(&m->tStart);

  return(m);
}

/* Recompute total number of splice sites and exons of current sequence */
void updateTotals(account *m,
                  packSites* allSites,
                  packSites* allSites_r,
                  packExons* allExons,
                  packExons* allExons_r)
{
  m->starts += allSites->nStartCodons;
  m->stops += allSites->nStopCodons;
  m->acc += allSites->nAcceptorSites;
  m->don += allSites->nDonorSites;

  m->starts_r += allSites_r->nStartCodons;
  m->stops_r += allSites_r->nStopCodons;
  m->acc_r += allSites_r->nAcceptorSites;
  m->don_r += allSites_r->nDonorSites;
   
  m->first += allExons->nInitialExons;
  m->internal += allExons->nInternalExons;
  m->terminal += allExons->nTerminalExons;
  m->single += allExons->nSingles;

  m->first_r += allExons_r->nInitialExons;
  m->internal_r += allExons_r->nInternalExons;
  m->terminal_r += allExons_r->nTerminalExons;  
  m->single_r += allExons_r->nSingles;  

  m->totalExons = 
    m->first + 
    m->first_r +
    m->internal + 
    m->internal_r +
    m->terminal +
    m->terminal_r +
    m->single +
    m->single_r;
}

/* Reset counters to the next DNA sequence */
void cleanAcc(account* m)
{
  
  /* ReInitialize */
  m->totalExons = 0;
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
  m->first_r = 0;
  m->internal_r = 0;
  m->terminal_r = 0;
  m->single_r = 0;
}

