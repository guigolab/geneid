/*************************************************************************
*                                                                        *
*   Module: BuildSort                                                    *
*                                                                        *
*   Sorting exons by donor position according to every class properties  *
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

/*  $Id: BuildSort.c,v 1.2 2001-12-18 15:18:32 eblanco Exp $  */

#include "geneid.h"

void BuildSort(dict *D,
               int nc[],
               int ne[],
               int UC[][MAXENTRY],
               int DE[][MAXENTRY],
               int nclass,
               long km[],
               exonGFF* **d,
               exonGFF *E,
               long nexons)
{
  long i,k;
  int j;
  int type;
  int class;
  char aux[MAXTYPE];
  
  /* Every exon will be classified into some sorting function (d) */
  /* Input exons are sorted by acceptor (left position) */
  for(i=0; i < nexons; i++)
    {
      aux[0]='\0';
      strcpy (aux, (E+i)->Type);
      strcat (aux, &((E+i)->Strand));
	  
      /* What's the type of exon? "Type+Strand" */
      type = getkeyDict(D,aux);
      
      /* Checking and getting exon type (dictionary) */
      if (type != NOTFOUND)
		{
		  /* Exon may belong to some upstream compatible classes (UC) */
		  for(j=0; j < nc[type]; j++)
			{
			  class = UC[type][j];
			  k = km[class]-1;
			  
			  /* Screening the exons sorted before: sorting by insertion */
			  while (k>=0 && (((E+i)->Donor->Position + (E+i)->offset2) 
							  < 
							  (d[class][k]->Donor->Position 
							   + d[class][k]->offset2)))  
				{
				  /* Shifting down previous exons */
				  d[class][k+1] = d[class][k];
				  k--;
				}
			  /* Insert new exon before the previously shifted exons */
			  d[class][k+1] = (E+i);
			  km[class]++;
			}
		} /* end if type found */
    } /* end forall exons */
}




