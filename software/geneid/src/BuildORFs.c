/*************************************************************************
*                                                                        *
*   Module: BuildSingles                                                 *
*                                                                        *
*   Using lists of splice sites, it builds single genes.                 *
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

/*  $Id: BuildORFs.c,v 1.1 2000-12-18 09:37:55 eblanco Exp $  */

#include "geneid.h"

extern long NUMEXONS;

long BuildSingles(site *Start, long nStarts, 
		  site *Stop, long nStops,
		  long cutPoint,
		  exonGFF *Exon) 
{
  int Frame;
  long i, j, js;
  long nSingles;
  
  /* main loop, for each Start codon */
  for (i=0, j=0, nSingles=0;
      (i < nStarts) && (j<nStops) &&
      (nSingles<(int)(NUMEXONS/RSINGL));
      i++)
    {
      Frame = (Start+i)->Position % 3;
   
      /* advance Stops to Start */
      while ( (j < nStops) && 
	      (((Stop+j)->Position+1) < (Start+i)->Position))
	j++;

      js=j;

      /* Only use some stops (stops before cutPoint) */
      /* find first Stop in Frame with Start: Stop closes Start frame */
      while ((js < nStops) && (((Stop+js)->Position+1) % 3 != Frame))
	js++;

      /* Only use some stops (stops before cutPoint) */
      /* first Stop after Start defines a single gen */
      if (js < nStops && (Stop+js)->Position >= cutPoint)
	{
	  /* Length rule about Single Genes */
	  if ( ((Stop+js)->Position + LENGTHCODON - (Start+i)->Position + 1)
	       >= 
	       SINGLEGENELENGTH)
	    {
	      (Exon + nSingles)->Acceptor = (Start+i);
	      (Exon + nSingles)->Donor = (Stop+js);
	      (Exon + nSingles)->Frame = 0;
	      (Exon + nSingles)->Remainder = 0;
	      strcpy((Exon + nSingles)->Type,"Single");
	      (Exon + nSingles)->Group = NOGROUP;
	      (Exon + nSingles)->evidence = 0;

	      nSingles++;
	    }
	}
    }

  if (nSingles>=(int)(NUMEXONS/RSINGL))
    printError("Too many predicted exons: Change RSINGL parameter");
  
  return(nSingles);
}
