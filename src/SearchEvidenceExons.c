/*************************************************************************
*                                                                        *
*   Module: SearchEvidenceExons                                          *
*                                                                        *
*   Searching evidence exons between (l1, l2).                           *
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

/*  $Id: SearchEvidenceExons.c,v 1.2 2001-02-07 17:37:50 eblanco Exp $  */

#include "geneid.h"

/* Looking for evidence exons to join with current split predictions */
long SearchEvidenceExons(packEvidence* pv, long l2)
{
  long i;
  
  /* Evidences with acceptor higher than l2 will be always ignored */
  i = pv->i1vExons;
  while ((i < pv->nvExons) && ((pv->vExons+i)->Acceptor->Position) < l2)
    i++;
    
  pv->i2vExons = i;

  return(pv->i2vExons - pv->i1vExons);
}

void SwitchCounters(packEvidence* pv)
{
  pv->i1vExons = pv->i2vExons;
}
