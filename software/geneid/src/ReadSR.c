/*************************************************************************
*                                                                        *
*   Module: ReadSR                                                       *
*                                                                        *
*   It reads similarity to protein regions in GFF format.                *
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

/*  $Id: ReadSR.c,v 1.2 2000-08-08 14:20:29 eblanco Exp $  */

#include "geneid.h"

void ReverseSR(packSR* sr)
{
  int x;
  int i, j;
  long p1Aux, p2Aux;
  double scoreAux;

  for(x=0; x < FRAMES; x++)
    for (i = 0, j = (sr->nRegions[FRAMES + x])-1; i <= j; i++, j--)
      {

	p1Aux = sr->sRegions[FRAMES + x][i].Pos1;
	p2Aux = sr->sRegions[FRAMES + x][i].Pos2;
	scoreAux = sr->sRegions[FRAMES + x][i].Score;

	sr->sRegions[FRAMES + x][i].Pos1 = sr->sRegions[FRAMES + x][j].Pos1;
	sr->sRegions[FRAMES + x][i].Pos2 = sr->sRegions[FRAMES + x][j].Pos2;
	sr->sRegions[FRAMES + x][i].Score = sr->sRegions[FRAMES + x][j].Score;

	sr->sRegions[FRAMES + x][j].Pos1 = p1Aux;
	sr->sRegions[FRAMES + x][j].Pos2 = p2Aux;
	sr->sRegions[FRAMES + x][j].Score = scoreAux;

      }
}

long ReadSR (char *FileName, packSR* sr, long LengthSequence)
{
  int i;
  FILE *file;
  char line[MAXLINE];
  char mess[MAXSTRING];
  long pos1, pos2;
  double score;
  char strand;
  short frame;

  if ((file=fopen(FileName, "r"))==NULL)
    printError("The similarity regions file cannot be open for read");
  
  /* Reset counters */
  for(i=0;i<STRANDS*FRAMES;i++)
    {
      sr->nRegions[i] = 0;
      sr->iRegions[i] = 0;
    }            

  /* Coments: line begins with # */
  /* gff format = "Name  Source  Type  Begin  End  Score  Strand  Frame */
  i = 0;
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      if(line[0]=='#')
	{
	  /* Skip this line */
	}
      else
	{
	  if ((sscanf(line, "%*s %*s %*s %ld %ld %lf %c %hd %*d",
		      &pos1,
		      &pos2,
		      &score,
		      &strand,
		      &frame) != 5) || (frame<1) || (frame>3))
	    {
	      sprintf(mess, "Error similarity regions: line %d\n",i);
	      printError(mess);
	    }
	  else
	    {
	      /* Allocating similarity region in packSR */
	      if (strand == '+')
		{
		  frame = (frame % 3);

		  sr->sRegions[frame][sr->nRegions[frame]].Pos1 = pos1;
		  sr->sRegions[frame][sr->nRegions[frame]].Pos2 = pos2;
		  sr->sRegions[frame][sr->nRegions[frame]].Score = score;
		  sr->nRegions[frame]++;
		}
	      else
		{
		  frame = (frame % 3);
		  sr->sRegions[FRAMES + frame][sr->nRegions[FRAMES + frame]].Pos1 = 
		    LengthSequence - pos2 + 1;
		  sr->sRegions[FRAMES + frame][sr->nRegions[FRAMES + frame]].Pos2 = 
		    LengthSequence - pos1 + 1;
		  sr->sRegions[FRAMES + frame][sr->nRegions[FRAMES + frame]].Score = score;
		  sr->nRegions[FRAMES + frame]++;
		}
	      sr->nTotalRegions++;
	    }
	}
      i++;
    }
  fclose(file);
  
  /* RVS strand SR must be reversed */
  ReverseSR(sr);

  return(sr->nTotalRegions);
}
