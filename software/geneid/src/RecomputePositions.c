/*************************************************************************
*                                                                        *
*   Module: RecomputePositions                                           *
*                                                                        *
*   Predicted signals in RVS: translation into FWD sense positions       *
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

/*  $Id: RecomputePositions.c,v 1.4 2006-05-25 14:23:58 talioto Exp $  */

#include "geneid.h"

/* Normalising positions on the reverse (and complemented) sequence */
void RecomputePositions(packSites* allSites, long l)
{
  long i;
  
  for (i=0; i < allSites->nStartCodons; i++)
    (allSites->StartCodons+i)->Position = 
      l-(allSites->StartCodons+i)->Position-1;

  for (i=0; i < allSites->nAcceptorSites; i++)
    (allSites->AcceptorSites+i)->Position = 
      l-(allSites->AcceptorSites+i)->Position-1;

  for (i=0; i < allSites->nU12gtagAcceptorSites; i++)
    (allSites->U12gtagAcceptorSites+i)->Position = 
      l-(allSites->U12gtagAcceptorSites+i)->Position-1;

  for (i=0; i < allSites->nU12gtagDonorSites; i++)
    (allSites->U12gtagDonorSites+i)->Position = 
      l-(allSites->U12gtagDonorSites+i)->Position-1;

  for (i=0; i < allSites->nU12atacAcceptorSites; i++)
    (allSites->U12atacAcceptorSites+i)->Position = 
      l-(allSites->U12atacAcceptorSites+i)->Position-1;

  for (i=0; i < allSites->nU12atacDonorSites; i++)
    (allSites->U12atacDonorSites+i)->Position = 
      l-(allSites->U12atacDonorSites+i)->Position-1;

  for (i=0; i < allSites->nDonorSites; i++)
    (allSites->DonorSites+i)->Position = 
      l-(allSites->DonorSites+i)->Position-1;

  for (i=0; i < allSites->nStopCodons; i++)
    (allSites->StopCodons+i)->Position = 
      l-(allSites->StopCodons+i)->Position-1;
}

   
