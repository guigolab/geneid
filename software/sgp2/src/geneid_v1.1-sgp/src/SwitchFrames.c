/*************************************************************************
*                                                                        *
*   Module: SwitchFrames                                                 *
*                                                                        *
*   Exchange frame and remainder from reverse exons                      *
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

/*  $Id: SwitchFrames.c,v 1.1 2003-09-10 14:53:34 gparra Exp $  */

#include "geneid.h"

/* Exchange frame and remainder in the input exons */
void SwitchFrames(exonGFF* e, long n)
{
  long i;
  int f;

  /* Exchange frame/rmd in reverse exons and reset the selected flag */
  for (i=0; i<n; i++)
	{      
	  if ((e+i)->Strand == '-')
		{
		  f=(e+i)->Frame;
		  (e+i)->Frame=(e+i)->Remainder;
		  (e+i)->Remainder=f;
		}

	  /* Mark exon as prediction in the current fragment */
	  (e+i)->selected = 0;
	}
}

/* Exchange frame and remainder in the sorted-by-donor exons only once */
/* Right now, d-array only contain exons from last fragment processing */
void SwitchFramesDa(packGenes* pg, int nclass)
{
  long i;
  long j;
  int f;

  /* Screening every class looking for exons... */
  for (i=0; i < nclass; i++) 
	{
      /* Traversing the list of exons in this class */
      for (j=0; j < pg->km[i]; j++)
		if (pg->d[i][j]->Strand == '-')
		  {
            /* Exchange frame/rmd only once */
            /* One exon might be in more than one list */
            if (!pg->d[i][j]->selected)
			  {
				f= pg->d[i][j]->Frame;
				pg->d[i][j]->Frame = pg->d[i][j]->Remainder;
				pg->d[i][j]->Remainder = f;

				/* Mark exon */
				pg->d[i][j]->selected = 1;
			  }
		  }
	}
}   

/* Restore original frame and remainder in exons from last fragment */
void SwitchFramesDb(packGenes* pg, int nclass)
{
  long i;
  long j;
  int f;

  /* Screening every class looking for exons... */
  for (i=0; i < nclass; i++)
    {
	  /* Traversing the list of exons in this class */
	  for (j=0; j < pg->km[i]; j++)
		if (pg->d[i][j]->Strand == '-')
		  {
			/* Exchange frame/rmd only once */
			/* Only exons from last fragment will have selected = 1 */
			if (pg->d[i][j]->selected)
			  {
				f = pg->d[i][j]->Frame;
				pg->d[i][j]->Frame = pg->d[i][j]->Remainder;
				pg->d[i][j]->Remainder = f;

				/* Mark exon */
				pg->d[i][j]->selected = 0;
			  }   
		  }
	}
}

/* Restore frame/remainder in reverse exons read from gff file */
void UndoFrames(exonGFF* e, long n)
{
  long i;
  int f;
  
  /* Undo frame/rmd change*/
  for (i=0;i<n;i++)  
    if ((e+i)->Strand == '-')
      {
		f=(e+i)->Frame;
		(e+i)->Frame=(e+i)->Remainder;
		(e+i)->Remainder=f;
      }
}
