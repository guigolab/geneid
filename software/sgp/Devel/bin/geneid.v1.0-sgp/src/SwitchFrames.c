/*************************************************************************
*                                                                        *
*   Module: SwitchFrames                                                 *
*                                                                        *
*   GenAmic needs exchanged frame/remainder from reverse exons.          *
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

/*  $Id: SwitchFrames.c,v 1.1 2000-07-05 08:28:26 eblanco Exp $  */

#include "geneid.h"

void SwitchFrames(exonGFF *e, long n, int when) 
{
  long i;
  int f;

  switch(when)
    {
    case 0:
      for (i=0;i<n;i++) 
	{
	  if ((e+i)->evidence)
            (e+i)->selected = 1; 
	  
	  /* Exchange frame/rmd only once */
	  if (((e+i)->Strand == '-') && !(e+i)->selected && !(e+i)->evidence)
	    {
	      f=(e+i)->Frame;
	      (e+i)->Frame=(e+i)->Remainder;
	      (e+i)->Remainder=f;
	      (e+i)->selected = 1;
	    }
	}
      break;
    case 1:
      for (i=0;i<n;i++) 
	{
	  /* Exchange frame/rmd only once */
	  if (((e+i)->Strand == '-') && (e+i)->selected)
	    {
	      f=(e+i)->Frame;
	      (e+i)->Frame=(e+i)->Remainder;
	      (e+i)->Remainder=f;
	      (e+i)->selected = 0;
	    }
	}
      break;
    }
}

void SwitchFramesDa(packGenes* pg, int nclass)
{
  long i;
  long j;
  int f;

  for (i=0; i < nclass; i++) 
  {
      for (j=0; j < pg->km[i]; j++)
         if (pg->d[i][j]->Strand == '-')
         {
	   /* Exchange frame/rmd only once */
	   if (!pg->d[i][j]->selected)
            {
	      f= pg->d[i][j]->Frame;
	      pg->d[i][j]->Frame = pg->d[i][j]->Remainder;
	      pg->d[i][j]->Remainder = f;
	      pg->d[i][j]->selected = 1;
	    }
	 }
   }
}   

void SwitchFramesDb(packGenes* pg, int nclass)
{
  long i;
  long j;
  int f;

  for (i=0; i < nclass; i++) 
    for (j=0; j < pg->km[i]; j++)       
      if (pg->d[i][j]->Strand == '-')
	{
	  /* Exchange frame/rmd only once */
	  if (pg->d[i][j]->selected)
	    {
	      f= pg->d[i][j]->Frame;
	      pg->d[i][j]->Frame = pg->d[i][j]->Remainder;
	      pg->d[i][j]->Remainder = f;
	      pg->d[i][j]->selected = 0;
	    }   
	}
}

void UndoFrames(exonGFF *e, long n)
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
