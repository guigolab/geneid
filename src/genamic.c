/*************************************************************************
*                                                                        *
*   Module: genamic                                                      *
*                                                                        *
*   Construction of genes using the algorithm GenAmic.                   *
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

/*  $Id: genamic.c,v 1.3 2000-07-12 12:09:20 eblanco Exp $  */

#include "geneid.h"

extern int GENEID;

void genamic(exonGFF *E, long nExons, packGenes* pg, gparam* gp)
{
  long i,j,j2;
  short frame,remainder;
  int h;
  char saux[MAXTYPE+1];
  int type,etype;
  long MaxDist;
  long MinDist;
  char mess[MAXSTRING];
  
  /* Starting process ... */
  printMess("Running Genamic ...");

  /* geneid sends this set of exonsGFF */
  sprintf(mess,"%ld exons GFF sent by geneid\n",nExons);
  printMess(mess);

  /* Frame and remainder of reverse-strand exons must be exchanged */
  if (GENEID)
    {
      SwitchFrames(E,nExons,0);
      SwitchFramesDa(pg, gp->nclass);
    }
 
  /* 1. Create a set of sorted by donor array of pointer to exons */
  printMess("Building Sorting Functions...");

  /* Build a set of sorting exons by donor functions */
  BuildSort(gp->D, gp->nc, gp->ne, 
	    gp->UC, gp->DE, gp->nclass, 
	    pg->km, pg->d, E, nExons);

  /* 2. Genamic Algorithm in Linear Time */
  printMess("Assembling Genes...");

  /* For each exon verify non-blocking and distances well defined */
  /* and after this, it could be assembly with the best gene before it */
  for(i=0; i< nExons; i++)
    {
      /* What is the type of this exon? */
      saux[0]='\0';
      strcpy (saux, (E+i)->Type);
      strcat (saux, &((E+i)->Strand));
      type = getkeyDict(gp->D,saux);
      frame = (E+i)->Frame;      
      
      (E+i)->PreviousExon = pg->Ghost;
      (E+i)->GeneScore = (E+i)->Score;
      
      if (type != NOTFOUND)
	{
	  /* For each class equivalent it builds the best gene 
	     ending in this exon if this class(etype) is not blocked */
	  for(h=0; h < gp->ne[type]; h++)
	    {
	      etype = gp->DE[type][h];
	      j = pg->je[etype];
	      MaxDist = gp->Md[etype];
	      MinDist = gp->md[etype];
	      
	      if ((pg->Ga[etype][frame]->Strand !='*') &&
		  (pg->Ga[etype][frame]->Donor->Position 
		   + 
		   pg->Ga[etype][frame]->offset2) 
		  < 
		  ((E+i)->Acceptor->Position + (E+i)->offset1) - MaxDist) 
		{
		  /* loop back searching the new best gene in this distance */
		  pg->Ga[etype][frame] = pg->Ghost; 
		  j2 = j-1;
		  while (j2>=0 && j2 < pg->km[etype]  && 
			 ((pg->d[etype][j2]->Donor->Position 
			   + pg->d[etype][j2]->offset2)
			  >= 
			  ((E+i)->Acceptor->Position + (E+i)->offset1)
			  - MaxDist))
		    {
		      remainder = pg->d[etype][j2]->Remainder;
		      if (pg->d[etype][j2]->GeneScore > 
			  pg->Ga[etype][remainder] -> GeneScore)
			pg->Ga[etype][remainder] = pg->d[etype][j2];
		      j2--;
		    }
		}

	      /* Loop forward: One scan over each donor-sort array */
	      while(j < pg->km[etype] && 
		    ((pg->d[etype][j]->Donor->Position 
		      + pg->d[etype][j]->offset2)
		     <= 
		     ((E+i)->Acceptor->Position + (E+i)->offset1) - MinDist))
		{
		  remainder = pg->d[etype][j]->Remainder;
		  if ((frame == remainder && 
		       ((pg->d[etype][j]->Donor->Position 
			 + pg->d[etype][j]->offset2)
			< 
			((E+i)->Acceptor->Position + (E+i)->offset1) - MaxDist)))
		    {
		      /* Skip this exon because max distance not ok */
		    }
		  else
		    {
		      if (pg->d[etype][j]->GeneScore > 
			  pg->Ga[etype][remainder]->GeneScore)
			pg->Ga[etype][remainder] = pg->d[etype][j];
		    }
		  j++;
		}
	      pg->je[etype] = j;
	      
	      /* Assembling the exon with the best compatible gene before */
	      /* Verify group rules if there are evidence exons */
	      if ((pg->Ga[etype][frame]->Group == (E+i)->Group ||
		   gp->block[etype] == NONBLOCK)
		  &&
		  ((pg->Ga[etype][frame]->GeneScore + (E+i)->Score) > (E+i)->GeneScore))
		{
		  (E+i)->GeneScore = pg->Ga[etype][frame]->GeneScore + (E+i)->Score;  
		  (E+i)->PreviousExon = pg->Ga[etype][frame];
		}
	    }
	  /* Updating the best gene */
	  if (((E+i)->GeneScore) > (pg->GOptim -> GeneScore))
	    pg->GOptim = (E+i);
	}
    }

  /* Undo the change between frame and remainder */
  if (GENEID)
    {
      SwitchFrames(E,nExons,1); 
      SwitchFramesDb(pg, gp->nclass);
    }
  else
    UndoFrames(E,nExons);
}

