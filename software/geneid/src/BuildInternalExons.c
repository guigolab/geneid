/*************************************************************************
*                                                                        *
*   Module: BuildInternalExons                                           *
*                                                                        *
*   Using lists of splice sites, it builds internal exons.               *
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

/*  $Id: BuildInternalExons.c,v 1.1 2000-07-07 08:21:05 eblanco Exp $  */

#include "geneid.h"

long BuildInternalExons(site *Acceptor, long nAcceptors, 
		        site *Donor, long nDonors,
		        site *Stop, long nStops,
		        int MaxDonors,
			exonGFF* Exon)
{
  struct iexon 
  {
    site *Acceptor;
    site *Donor;
    int Frame[3];
  } *LocalExon;
  
  int nLocalExons,
    LowestLocalExon;
  double LowestLocalScore;
  
  int Frame[3];

  long i,
       j=0, js,
       k=0, ks;

  int f,
      l,
      ll;

  long nExon=0;
  LocalExon = (struct iexon *) calloc(MaxDonors, sizeof(struct iexon));
  
  /* main loop, for each Acceptor Site */
  for (i=0; (i < nAcceptors) && (nExon<NUMEXONS); i++) 
    {
      /* Open the windows... */
      for (f=0;f<3;f++) 
	Frame[f]=1;
      
      nLocalExons = 0;
      LowestLocalScore = DBL_MAX;
      
      /* advance Stop to Acceptor+i */
      while (((Stop+j)->Position+1 < (Acceptor+i)->Position) && (j < nStops))
	j++;
      
      js=j;
      
      /* advance Donor to Acceptor+i */
      while (((Donor+k)->Position < (Acceptor+i)->Position + EXONLENGTH) 
	     && (k < nDonors))
	k++;
      ks=k;
      
      /* find Stop in first Frame still open */
      while ((Frame[0]==1 || Frame[1]==1 || Frame[2]==1) && (ks < nDonors)) 
	{
	  /* if frame is open... */
	  if (js == nStops ||
	      (Frame[f=((Stop+js)->Position - (Acceptor+i)->Position + 1) % 3])) 
	    {
	      /* ...every Donor between Acceptor+i and Stop+js defines an exon */
	      /* keep  only the top scoring MaxDonors */
	      while ((js == nStops || (Donor+ks)->Position < (Stop+js)->Position + 1 + 2)
		     && (ks < nDonors))
		{
		  if (nLocalExons < MaxDonors) 
		    {
		      (LocalExon+nLocalExons)->Acceptor=(Acceptor+i);
		      (LocalExon+nLocalExons)->Donor=(Donor+ks);
		      for (ll=0;ll<3;ll++)
			(LocalExon+nLocalExons)->Frame[ll]=Frame[ll];
		      
		      if ((Donor+ks)->Score < LowestLocalScore) 
			{
			  LowestLocalScore = (Donor+ks)->Score;
			  LowestLocalExon = nLocalExons;
			}
		      nLocalExons++;
		    }
		  else 
		    {
		      if ((Donor+ks)->Score >= LowestLocalScore) 
			{
			  for (l=LowestLocalExon;l<nLocalExons-1;l++)
			    LocalExon[l]=LocalExon[l+1];
			  (LocalExon+nLocalExons-1)->Acceptor=(Acceptor+i);
			(LocalExon+nLocalExons-1)->Donor=(Donor+ks);
			for (ll=0;ll<3;ll++)
			  (LocalExon+nLocalExons-1)->Frame[ll]=Frame[ll];
			
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
	      Frame[f]=0;
	    }
	  (js<nStops)?
	    js++:js;
	}

      /* Save predicted exons with current Acceptor */
      for (l=0;(l<nLocalExons) && (nExon<NUMEXONS);l++) 
	{
	  for (ll=0;(ll<3) && (nExon<NUMEXONS);ll++)
	    if ((LocalExon+l)->Frame[ll]) {
	      (Exon+nExon)->Acceptor = (LocalExon+l)->Acceptor;
	      (Exon+nExon)->Donor = (LocalExon+l)->Donor;
	      (Exon+nExon)->Frame = ll;
	      (Exon+nExon)->Remainder = ((Exon+nExon)->Donor->Position -
					 ((Exon+nExon)->Acceptor->Position +
					  (Exon+nExon)->Frame)+ 1 ) % 3;
	      (Exon+nExon)->Remainder = (3 - (Exon+nExon)->Remainder) % 3;
	      strcpy((Exon+nExon)->Type,"Internal");
	      (Exon+nExon)->Group = NOGROUP;
	      (Exon+nExon)->evidence = 0;
	      nExon++;
	    }
	}
    }
  
  if (nExon==NUMEXONS)
    printError("Too many predicted exons: Change NUMEXONS parameter");
  
  free(LocalExon);
  return(nExon);
}


    









