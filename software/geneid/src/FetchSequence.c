/*************************************************************************
*                                                                        *
*   Module: FetchSequence                                                *
*                                                                        *
*   Put every base to upper letter and reverse sequence                  *
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

/*  $Id: FetchSequence.c,v 1.1 2000-07-05 08:15:23 eblanco Exp $  */

#include "geneid.h"

int complement(int c)
{

  switch (c) {
      case 'A':
        return('T');
      case 'C':
        return('G');
      case 'G':
        return('C');
      case 'T':
        return('A');
      default:
        return('N');
      }
}

/* Puts upper letters and reverse/complement a DNA sequence */
long FetchSequence(char *s, char* r)
{
  long i,j;
  long l;

  l = strlen(s);
  for (i = 0, j = l-1; i <= j; i++, j--)
    {    
      /* Upper case */
      if (s[i] > 96)
         s[i] = s[i] - 32;
      if (s[j] > 96)
         s[j] = s[j] - 32;

      /* Reverse */
      r[i]=complement(s[j]);
      r[j]=complement(s[i]);
    }
  r[l] = '\0';

 return (l);
}

/* Reverse a split of the original sequence */
void ReverseSubSequence(long p1, long p2, char* s, char* r)
{
  long i,j;
  long l;

  l = p2 - p1 + 1;
  for (i = p1, j = p1 + (l-1); i <= j; i++, j--)
    {    
      /* Reverse */
      r[i-p1]=complement(s[j]);
      r[j-p1]=complement(s[i]);
    }
}


