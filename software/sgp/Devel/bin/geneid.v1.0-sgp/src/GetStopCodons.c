/*************************************************************************
*                                                                        *
*   Module: GetStopCodons                                                *
*                                                                        *
*   Stop codons prediction using Position Frequency Matrices.            *
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

/*  $Id: GetStopCodons.c,v 1.2 2000-10-03 10:52:13 jabril Exp $  */

#include "geneid.h"

extern int TRANS[];
extern long NUMSITES;

long GetStopCodons(char *s, profile *p, 
		   site *sc, 
		   long l1, 
		   long l2) 
{
  long ns=0,
       is=0;
  double score;
  int i;
  
  static char *stop[] = 
  {
    "TAA",
    "TAG",
    "TGA"
  };

  char codon[LENGTHCODON];
  long left, 
       right;
  int index;

  /* left and right are the true boundaries of prediction */
  left  = MAX(0+p->order,l1 - p->offset);
  right = l2 - p->offset;
  s += left;

  if (p->order==0)
    {
      /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right - left + 1) && (ns<NUMSITES))  
	{
	  strncpy(codon,s+p->offset+1,LENGTHCODON); 
	  codon[LENGTHCODON]='\0';
	  
	  if (!strcmp(codon,stop[0]) ||
	      !strcmp(codon,stop[1]) ||
	      !strcmp(codon,stop[2])) 
	    {
	      score=0.0;
	      for (i=0;i<p->dimension;i++) 
		{
		  /* i is the position inside the region */
		  index = TRANS[(int)(*(s + i))];
		  if (index >= p->dimensionTrans)
		    score = score + -INFI;
		  else
		    score = score + p->transitionValues[i][index];
		}
	      
	      if (score >= p->cutoff) 
		{
		  /* Position given is last coding position before Stop */ 
		  sc[ns].Position=left + is + p->offset;
		  sc[ns].Score=score;
		  ns++;
		}
	    } 
	  is++;
	  s++;
	}
    }   
  else if (p->order==1)
    {
      /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right - left + 1) && (ns<NUMSITES))  
	{
	  strncpy(codon,s+p->offset+1,LENGTHCODON); 
	  codon[LENGTHCODON]='\0';
	  
	  if (!strcmp(codon,stop[0]) ||
	      !strcmp(codon,stop[1]) ||
	      !strcmp(codon,stop[2])) 
	    {
	      score=0.0;
	      for (i=0;i<p->dimension;i++) 
		{
		  /* i is the position inside the region */
		  index = 4*TRANS[(int)(*(s + i -1))] + TRANS[(int)(*(s + i))];
		  if (index >= p->dimensionTrans)
		    score = score + -INFI;
		  else
		    score = score + p->transitionValues[i][index];
		}
	      
	      if (score >= p->cutoff) 
		{
		  /* Position given is last coding position before Stop */ 
		  sc[ns].Position=left + is + p->offset;
		  sc[ns].Score=score;
		  ns++;
		}
	    } 
	  is++;
	  s++;
	}     
    }
  else
    {
     /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right - left + 1) && (ns<NUMSITES))  
	{
	  strncpy(codon,s+p->offset+1,LENGTHCODON); 
	  codon[LENGTHCODON]='\0';
	  
	  if (!strcmp(codon,stop[0]) ||
	      !strcmp(codon,stop[1]) ||
	      !strcmp(codon,stop[2])) 
	    {
	      score=0.0;
	      for (i=0;i<p->dimension;i++) 
		{
		  /* i is the position inside the region */
		  index = OligoToInt(s + i - p->order , p->order+1); 
		  if (index >= p->dimensionTrans)
		    score = score + -INFI;
		  else
		    score = score + p->transitionValues[i][index];
		}
	      
	      if (score >= p->cutoff) 
		{
		  /* Position given is last coding position before Stop */ 
		  sc[ns].Position=left + is + p->offset;
		  sc[ns].Score=score;
		  ns++;
		}
	    } 
	  is++;
	  s++;
	} 
    }

  /* Get remaining stops until the end of sequence. Set score to 0 for those */
  if (!(s+p->dimension))
    {
      s=(s-is);
      is+=p->offset;
      
      while (*(s+is) && (ns<NUMSITES))  
	{
	  strncpy(codon,(s+is),3); 
	  codon[3]='\0';
	  if (!strcmp(codon,stop[0]) || !strcmp(codon,stop[1]) || 
	      !strcmp(codon,stop[2])) 
	    {
	      /* Position given is last coding position before Stop */
	      sc[ns].Position=is-1;  
	      sc[ns].Score=0;
	      ns++;
	    } 
	  is++;
	}
    }
  
  if (ns==NUMSITES)
    printError("Too many predicted sites: Change NUMSITES parameter");

  return(ns);
}
