/*************************************************************************
*                                                                        *
*   Module: SetRatios                                                    *
*                                                                        *
*   Stablish maximum allowed number of signals and exons per split       *
*                                                                        *
*   This file is part of the geneid 1.4 distribution                     *
*                                                                        *
*     Copyright (C) 2006 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
*                          Tyler   ALIOTO                                * 
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

#include "geneid.h"

extern int BEG;

/* NUMSITES, NUMEXONS:
From the defined values RSITES, REXONS (geneid.h) and length of sequence,
an estimation for the amount of predicted signals and exons (any type) is
computed. Historically this sized the fixed site/exon arrays up front; now
that those arrays grow on demand (see GrowSiteArray/GrowExonArray), NUMSITES
and NUMEXONS are only a rough estimate left over for beggar.c's -B memory
report -- they no longer bound or size anything real. */

/* MAXBACKUPSITES, MAXBACKUPEXONS:
From the defined values RBSITES, RBEXONS (geneid.h) and length of sequence,
an estimation for the amount of signals and exons (any type), necessary
to restore the prediction between 2 splits, is computed. These two ARE still
real: MAXBACKUPEXONS sizes the dumpster's hash table (MAXDUMPHASH, see
DumpHash.c), and both being nonzero is the "multi-split, dumpster active"
flag checked in BackupGenes.c -- the dumpster's own backup arrays are
chunked and grow on demand (see BackupGenes.c), but the hash stays a fixed
table sized from this estimate. */

void SetRatios(long* NUMSITES,
               long* NUMEXONS,
               long* MAXBACKUPSITES,
               long* MAXBACKUPEXONS,
               long L)
{
  /* L is the estimated length of input DNA sequence */
  /* LENGTHSi is the length splitting */
  if (L < LENGTHSi)
    {
      /* Short sequences processed as a whole: only one split, so there is
	 no inter-split state to restore -- the 0s below are read as the
	 "dumpster inactive" flag in BackupGenes.c. */
      *NUMSITES = L / RSITES + BASEVALUESITES_SHORT;
      *NUMEXONS = L / REXONS + BASEVALUEEXONS_SHORT;

      /* There is no need to divide the sequence */
      *MAXBACKUPSITES = 0;
      *MAXBACKUPEXONS = 0;
    }
  else
    {
      /* Long sequences must be divided into several fragments */
      *NUMSITES = LENGTHSi / RSITES;
      *NUMEXONS = LENGTHSi / REXONS;

      /* Information inter-split predictions must be saved */
      *MAXBACKUPSITES = (L / RBSITES) + BASEVALUESITES_LARGE;
      *MAXBACKUPEXONS = (L / RBEXONS) + BASEVALUEEXONS_LARGE;
    }

}
