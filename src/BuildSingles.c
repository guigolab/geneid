/*************************************************************************
*                                                                        *
*   Module: BuildSingles                                                 *
*                                                                        *
*   From start and stop codons, to build single gene exons               *
*                                                                        *
*   This file is part of the geneid 1.4 distribution                     *
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

/*  $Id: BuildSingles.c,v 1.6 2011-01-13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic exons (divided by RSINGL) */

long BuildSingles(site *Start, long nStarts, 
                  site *Stop, long nStops,
                  long cutPoint, 
				  char* Sequence,
                  exonGFF **ExonP, long* capP)
{
  int Frame;
  long i, j, js;

  /* Growable output array (grows on demand; capacity tracked by caller) */
  exonGFF* Exon = *ExonP;
  long cap = *capP;

  /* Final number of predicted single genes exons */
  long nSingles;

  /* Main loop, for each Start codon searching the first Stop in frame */
  for (i=0, j=0, nSingles=0;
	   (i < nStarts) && (j<nStops);
	   i++)
    {
      Frame = (Start+i)->Position % 3;
	  
      /* Skip previous Stops to Start */
      while ( (j < nStops) && (((Stop+j)->Position+1) < (Start+i)->Position))
        j++;
	  
      /* Save counter j for the next iteration */
      js=j;
	  
      /* Skip Stops not in frame with the current Start */
      while ((js < nStops) && (((Stop+js)->Position+1) % 3 != Frame))
        js++;
	  
      /* CutPoint: to preserve sorted exons between fragments */
      if (js < nStops && (Stop+js)->Position >= cutPoint)
		{
		  /* LENGTH rule about Single Genes */
		  if ( ((Stop+js)->Position + LENGTHCODON - (Start+i)->Position + 1)
			   >= SINGLEGENELENGTH)
			{
			  GrowExonArray(&Exon,&cap,nSingles+1);
			  (Exon + nSingles)->Acceptor = (Start+i);
			  (Exon + nSingles)->Donor = (Stop+js);
			  (Exon + nSingles)->Frame = 0;
			  (Exon + nSingles)->Remainder = 0;
			  strcpy((Exon + nSingles)->Type,"Single");
			  strcpy((Exon + nSingles)->Group,NOGROUP);
			  (Exon + nSingles)->evidence = 0;
			  
			  nSingles++;
			}
		}
    }

  *ExonP = Exon;
  *capP = cap;
  return(nSingles);
}
