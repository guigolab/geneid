/*************************************************************************
*                                                                        *
*   Module: SearchEvidenceExons                                          *
*                                                                        *
*   Extracting evidence exons between (l1, l2) for the current split     *
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

/*  $Id: SearchEvidenceExons.c,v 1.1 2003-09-10 14:53:34 gparra Exp $  */

#include "geneid.h"

/* Looking for annotations to include with current ab initio predictions */
void SearchEvidenceExons(packEvidence* pv, long l2)
{
  long i;
  
  /* Evidences with acceptor higher than l2 will be always ignored */
  i = pv->i1vExons;
  while ((i < pv->nvExons) &&
		 ((pv->vExons+i)->Acceptor->Position + (pv->vExons+i)->offset1) <= l2)
    i++;
  
  pv->i2vExons = i;
  
  pv->ivExons = pv->i2vExons - pv->i1vExons;
}

/* Processing next block of annotations */
void SwitchCounters(packEvidence* pv)
{
  pv->i1vExons = pv->i2vExons;
}

/* Reset counters in packEvidence for the next input sequence */
void resetEvidenceCounters(packEvidence* pv)
{
  pv->i1vExons = 0;
  pv->i2vExons = 0;
}



