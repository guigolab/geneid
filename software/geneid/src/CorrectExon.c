/*************************************************************************
*                                                                        *
*   Module: CorrectExon                                                  *
*                                                                        *
*   Recompute splice sites positions according to its type               *
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

/*  $Id: CorrectExon.c,v 1.1 2000-07-05 08:10:34 eblanco Exp $  */

#include "geneid.h"

void CorrectExon(exonGFF *e) 
{
  /* Correct positions of Begin & End of this Exon */
  /* Stop Codon included into the exon */
  e->offset1 =
    ((!strcmp(e->Type,"Terminal") || !strcmp(e->Type,"Single"))
	     && e->Strand =='-') ? 
    -LENGTHCODON + COFFSET : +COFFSET;
  
  e->offset2 = 
    ((!strcmp(e->Type,"Terminal") || !strcmp(e->Type,"Single"))
	     && e->Strand =='+')? 
    LENGTHCODON + COFFSET : +COFFSET;
}



