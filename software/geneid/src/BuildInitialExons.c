/*************************************************************************
*                                                                        *
*   Module: BuildInitialExons                                            *
*                                                                        *
*   Using lists of splice sites, it builds initial exons                 *
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

/*  $Id: BuildInitialExons.c,v 1.2 2000-08-08 14:15:10 eblanco Exp $  */

#include "geneid.h"

extern long NUMEXONS;

long BuildInitialExons(site *Start, long nStarts, 
		       site *Donor, long nDonors,
		       site *Stop, long nStops,
		       int MaxDonors,
		       exonGFF *Exon) 
{
  exonGFF *LocalExon;
  int nLocalExons, LowestLocalExon;
  double LowestLocalScore;
  int Frame;
  long i, j=0, js, k=0, ks;
  int l;
  long nExon;
  
  LocalExon = (exonGFF*) calloc(MaxDonors, sizeof(exonGFF));
  
  /* main loop, for each Start codon */
  for (i = 0, nExon = 0;
       (nExon<(int)(NUMEXONS/RFIRST)) && (i < nStarts); 
       i++)
  { 
    Frame = (Start+i)->Position % 3;
    nLocalExons = 0;
    LowestLocalScore = DBL_MAX;
    
    /* advance Stops to Start */
    while (((Stop+j)->Position+1 < (Start+i)->Position) && (j < nStops))
      j++;
    js=j;
    
    /* find first Stop in Frame with Start */
    while ((((Stop+js)->Position+1) % 3 != Frame) && (js < nStops))
      js++;
    
    /* advance Donors to Start */
    while (((Donor+k)->Position < (Start+i)->Position+2) && (k < nDonors))
      k++;
    ks=k;
    
    /* every Donor between Start and Stop[js] defines an initial exon */
    /* keep only the top scoring MaxDonors */

    while ((js == nStops || (Donor+ks)->Position < (Stop+js)->Position + 1 + 2) 
	   && (ks < nDonors)) 
      {
	if (nLocalExons < MaxDonors) 
	  {
	    (LocalExon+nLocalExons)->Acceptor=(Start+i);
	    (LocalExon+nLocalExons)->Donor=(Donor+ks);
	    if ((Donor+ks)->Score < LowestLocalScore) 
	      {
		LowestLocalScore = (Donor+ks)->Score;
		LowestLocalExon = nLocalExons;
	      }
	    nLocalExons++;
	  }
	else 
	  {
	    if ((Donor+ks)->Score > LowestLocalScore) {
	      for (l=LowestLocalExon;l<nLocalExons-1;l++)
		LocalExon[l]=LocalExon[l+1];
	    (LocalExon+nLocalExons-1)->Acceptor=(Start+i);
	    (LocalExon+nLocalExons-1)->Donor=(Donor+ks);
	    
	    LowestLocalExon = 0;
	    LowestLocalScore = (LocalExon+0)->Donor->Score;
	    for (l=1;l<nLocalExons;l++)
	      if ((LocalExon+l)->Donor->Score < LowestLocalScore) 
		{
		  LowestLocalScore = (LocalExon+l)->Donor->Score;
		  LowestLocalExon = l;
		}
	    }
	  }
	ks++;
      }
    
    /* put LocalExon into InitialExons. Within Start sort by Donor */  
    for (l=0;(l<nLocalExons) && (nExon<(int)(NUMEXONS/RFIRST));l++) 
      {
	Exon[nExon] = LocalExon[l];
	(Exon+nExon)->Frame = 0;
	(Exon+nExon)->Remainder = ((Exon+nExon)->Donor->Position -
				   (Exon+nExon)->Acceptor->Position + 1 ) % 3;
	(Exon+nExon)->Remainder = (3 - (Exon+nExon)->Remainder) % 3; 
	strcpy((Exon+nExon)->Type,"First");
	(Exon+nExon)->Group = NOGROUP;
	(Exon+nExon)->evidence = 0;
	nExon++;
    }
  }

  if (nExon>=(int)(NUMEXONS/RFIRST))
   printError("Too many predicted exons: Change RFIRST parameter");

  free(LocalExon);
  return(nExon);
}
