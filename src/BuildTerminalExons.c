/*************************************************************************
*                                                                        *
*   Module: BuildTerminalExons                                           *
*                                                                        *
*   Using lists of splice sites, it builds terminal exons                *
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

/*  $Id: BuildTerminalExons.c,v 1.2 2000-08-08 14:16:16 eblanco Exp $  */

#include "geneid.h"

extern long NUMEXONS;

long BuildTerminalExons (site *Acceptor, long nAcceptors, 
		         site *Stop, long nStops,
			 long LengthSequence,
			 exonGFF* Exon,
			 long cutPoint)
{
  int Frame[3];
  long StartSearch;
  long iA=0,
       i, f,
       j=0, js;

  long nExon=0;

  /* advance Acceptor */
  StartSearch = 0;
  while((Acceptor+iA)->Position < StartSearch)
    iA++;
  
  while((Stop+j)->Position+1 < StartSearch)
    j++;
  
  /* main loop, for each Acceptor next Stop in every Frame defines an exon */
  for (i=iA;(i<nAcceptors) && (nExon<(int)(NUMEXONS/RTERMI));i++) 
    {
      /* Open the windows... */
      for (f=0;f<3;f++)
	Frame[f]=1;

      /* advance Stops to Acceptor */
      while (((Stop+j)->Position+1 < (Acceptor+i)->Position) && (j < nStops))
	j++;
      js=j;
      
      /* Find Stops in Frame still open */
      while ((Frame[0]==1 || Frame[1]==1 || Frame[2]==1)
	     && (js < nStops)
	     && (nExon<(int)(NUMEXONS/RTERMI)))
	{
	  if (Frame[f=((Stop+js)->Position - (Acceptor+i)->Position + 1) % 3])
	    {
	      if ((Stop+js)->Position >= (Acceptor+i)->Position &&
		  (Stop+js)->Position >= cutPoint)
		{
		  (Exon+nExon)->Acceptor=(Acceptor+i);
		  (Exon+nExon)->Donor=(Stop+js);
		  (Exon+nExon)->Frame = f;
		  (Exon+nExon)->Remainder = 0; 
		  strcpy((Exon+nExon)->Type,"Terminal");
		  (Exon+nExon)->Group = NOGROUP;
		  (Exon+nExon)->evidence = 0;
		  nExon++;
		}
	      Frame[f]=0;
	    } 
	  js++;
	}     
    }    

  if (nExon >=(int)(NUMEXONS/RTERMI))
    printError("Too many predicted exons: Change RTERMI parameter");
  
  return(nExon);
}
