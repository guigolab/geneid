/*************************************************************************
*                                                                        *
*   Module: genamic                                                      *
*                                                                        *
*   Assembling genes from the input set of exons                         *
*                                                                        *
*   This file is part of the geneid 1.1 Distribution                     *
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

/*  $Id: genamic.c,v 1.4 2001-12-18 16:22:03 eblanco Exp $  */

#include "geneid.h"

/* Complete gene prediction (sites and exons) or only assembling */
extern int GENEID;

void genamic(exonGFF* E, long nExons, packGenes* pg, gparam* gp)
{
  long i,j,j2;
  short frame,remainder;
  int h;
  char saux[MAXTYPE];
  int type,etype;
  long MaxDist;
  long MinDist;
  char mess[MAXSTRING];
  
  /* 0. Starting process ... */
  printMess("-- Running gene assembling (genamic) --");

  /* geneid sends this set of exonsGFF together with the Ga and d-array*/
  sprintf(mess,"%ld exons GFF received from geneid",nExons);
  printMess(mess);

  /* Frame/remainder of reverse-strand exons must be exchanged */
  if (GENEID)
    {
      /* Exons from the current fragment of sequence */
      SwitchFrames(E,nExons);
	  
      /* Exons from the last fragment of sequence */
      SwitchFramesDa(pg,gp->nclass);
    }
  
  /* 1. Create a set of sorted by donor array of pointer to exons */
  printMess("Sorting exons by donor...");
  
  /* Build a set of sorting exons by donor functions */
  BuildSort(gp->D, gp->nc, gp->ne, 
			gp->UC, gp->DE, gp->nclass, 
			pg->km, pg->d, E, nExons);
  
  /* 2. Genamic Algorithm in linear time (size of input) */
  printMess("Assembling Genes...");

  /* For every exon verify non-blocking and distances well defined */
  /* and after this, it might be assemble with the best gene computed */
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
		  /* For every equivalent class building the best gene ending with it */
		  for(h=0; h < gp->ne[type]; h++)
			{
			  etype = gp->DE[type][h];
			  j = pg->je[etype];
			  MaxDist = gp->Md[etype];
			  MinDist = gp->md[etype];
			  
			  /* Checking maximum distance allowed requirement */
			  if ((MaxDist != INFI) &&
				  (pg->Ga[etype][frame]->Strand !='*') &&
				  (pg->Ga[etype][frame]->Donor->Position 
				   + 
				   pg->Ga[etype][frame]->offset2) 
				  < 
				  ((E+i)->Acceptor->Position + (E+i)->offset1) - MaxDist) 
				{
				  /* loop backward searching another best gene matching MAX distance */
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
			  /* while minimum distance allowed requirement is OK */
			  /* Update best partial genes between current and previous exon */
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
					  /* Updating best partial assembled gene */
					  if (pg->d[etype][j]->GeneScore > pg->Ga[etype][remainder]->GeneScore)
						{
						  pg->Ga[etype][remainder] = pg->d[etype][j];
						}													
					}
				  j++;
				}
			  pg->je[etype] = j;
			  
			  /* Assembling the exon with the best compatible gene before it */
			  /* Verify group rules if there are evidence exons (annotations) */
			  if ((!(strcmp(pg->Ga[etype][frame]->Group,(E+i)->Group))
				   || gp->block[etype] == NONBLOCK)
				  &&
				  ((pg->Ga[etype][frame]->GeneScore + (E+i)->Score) > (E+i)->GeneScore))
				{
				  if (((E+i)->Strand == '+') && 
					  ((pg->Ga[etype][frame]->rValue == 1 && (E+i)->lValue == 1)
					   ||
					   (pg->Ga[etype][frame]->rValue == 2 && (E+i)->lValue == 2)
					   ||
					   (pg->Ga[etype][frame]->rValue == 3 && ((E+i)->lValue == 2 || (E+i)->lValue == 3))))
					{
					  /* FWD: Avoiding building a stop codon */
					}
				  else
					{
					  if (((E+i)->Strand == '-') && 
						  ((pg->Ga[etype][frame]->lValue == 1 && (E+i)->rValue == 1)
						   ||
						   (pg->Ga[etype][frame]->lValue == 2 && (E+i)->rValue == 2)
						   ||
						   ((pg->Ga[etype][frame]->lValue == 2 || pg->Ga[etype][frame]->lValue == 3) && (E+i)->rValue == 3)))
						{
						  /* RVS: Avoiding building a stop codon */
						}
					  else
						{
						  (E+i)->GeneScore = pg->Ga[etype][frame]->GeneScore + (E+i)->Score;  
						  (E+i)->PreviousExon = pg->Ga[etype][frame];
						}
					}
				}
			}

		  /* Updating the best gene assembled (final gene) */
		  if (((E+i)->GeneScore) > (pg->GOptim -> GeneScore))
			pg->GOptim = (E+i);
		}
    }

  /* 3. Undo the change between frame and remainder */
  if (GENEID)
    {
      /* Exons predicted in the current fragment of sequence */
      SwitchFrames(E,nExons);
	  
      /* Only changing exons predicted in the last fragment of sequence */
      SwitchFramesDb(pg,gp->nclass);
    }
  else
    UndoFrames(E,nExons);
  
  /* Finishing process */
  printMess("-- Finishing gene assembling (genamic) --");
}

