/*************************************************************************
*                                                                        *
*   Module: ReadSR                                                       *
*                                                                        *
*   Loading Similarity to protein Regions (SR) in GFF format             *
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

/*  $Id: ReadSR.c,v 1.3 2001-12-18 15:58:18 eblanco Exp $  */

#include "geneid.h"

/* Input SRs from external file which need not to be sorted */
long ReadSR (char* FileName, packSR* sr, long LengthSequence)
{
  int i;
  FILE* file;
  char line[MAXLINE];
  char lineCopy[MAXLINE];
  char *line1;
  char *line2;
  char *line3;
  char *line4;
  char *line5;
  char *line6;
  char *line7;
  char *line8;
  char *line9;
  
  /* If frame = '.' then make three copies of current SR (3 frames) */
  int three;
  
  char mess[MAXSTRING];
  long pos1, pos2;
  double score;
  char strand;
  short frame;
  char c;
  
  /* Coments: line begins with # */
  /* gff format = Name  Source  Type  Begin  End  Score  Strand  Frame */
  /* group is allowed but skipped */
  
  /* 0. Trying to open the SR file */
  if ((file=fopen(FileName, "r"))==NULL)
    printError("The similarity regions file can not be opened to read");
  
  /* 1. Read SRs */
  i = 0;
  three = 0;
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      /* 2.a. Comment or empty line: forget it */
      if (line[0]=='#' || line[0]=='\n')
		printMess("Skipping comment line in SRs file");
      else
		{
          /* Backup copy of the line to show error messages */
          strcpy(lineCopy,line);
		  
          /* Extracting GFF features */
          line1 = (char *) strtok(line,"\t");
          line2 = (char *) strtok(NULL,"\t");
          line3 = (char *) strtok(NULL,"\t"); 
          line4 = (char *) strtok(NULL,"\t");
          line5 = (char *) strtok(NULL,"\t");
          line6 = (char *) strtok(NULL,"\t");
          line7 = (char *) strtok(NULL,"\t");
          line8 = (char *) strtok(NULL,"\t");
          line9 = (char *) strtok(NULL,"\n");
		  
		  /* There are 8 mandatory columns and the last one is optional */
		  if (line1 == NULL || line2 == NULL || line3 == NULL ||
			  line4 == NULL || line5 == NULL || line6 == NULL ||
			  line7 == NULL || line8 == NULL)
			{
			  sprintf(mess,"Wrong GFF format in SRs (number of records):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 1/2/3/9. Sequence/Source/SR/group records not used anymore */
		  
		  /* 4. Starting position */
		  if (sscanf(line4,"%ld",&pos1) != 1)
			{
			  sprintf(mess,"Wrong GFF format in SRs (starting position):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 5. Finishing position */
		  if (sscanf(line5,"%ld",&pos2) != 1)
			{
			  sprintf(mess,"Wrong GFF format in SRs (finishing position):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 6. SR score */
		  if (sscanf(line6,"%lf",&score) != 1)
			{
			  sprintf(mess,"Wrong GFF format in SRs (score):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 7. Strand (reading sense) [+|-] */
		  if ((sscanf(line7,"%c",&strand)!= 1) ||
			  ((strand != '+') && (strand != '-')))
			{       
			  sprintf(mess,"Wrong GFF format in SRs (strand):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 8. Frame = integer or '.' */
		  if (sscanf(line8,"%hd",&frame) != 1)
			{
			  /* Is it a dot? */
			  if ((sscanf(line8,"%c",&c)!= 1) || (c!='.'))
				{
				  sprintf(mess,"Wrong GFF format in SRs (frame):\n-->%s\n",lineCopy);
				  printError(mess);
				}
			  /* make three copies for this SR */
			  three = 1;
			}
		  else
			{
			  /* Checking input frame between 1..3 */
			  if ((frame < 1) || (frame > 3))
				{
				  sprintf(mess,"Wrong GFF value in SRs ((blast) frame not between 1..3):\n-->%s\n",lineCopy);
				  printError(mess);
				}
			}
		  
		  /* Allocating similarity region in packSR */
		  /* a) Forward sense */
		  if (strand == '+')
			{
			  /* a1. Unknown frame: 3 copies */
			  if (three)
				{
				  /* Replication: 3 SRs for the same coordinates */
				  for (frame=0; frame < FRAMES; frame++)
					{
					  /* Setting SR data */
					  sr->sRegions[frame][sr->nRegions[frame]].Pos1 = pos1;
					  sr->sRegions[frame][sr->nRegions[frame]].Pos2 = pos2;
					  sr->sRegions[frame][sr->nRegions[frame]].Score = score;

					  sr->nRegions[frame]++;
					}
				}
			  /* a2. Known frame */
			  else
				{
				  /* Convert blast frame 3 into frame 0 to store it */
				  frame = (frame % 3);

				  /* Setting SR data */
				  sr->sRegions[frame][sr->nRegions[frame]].Pos1 = pos1;
				  sr->sRegions[frame][sr->nRegions[frame]].Pos2 = pos2;
				  sr->sRegions[frame][sr->nRegions[frame]].Score = score;

				  sr->nRegions[frame]++;
				}
			} /* end FWD */
		  /* b) Reverse sense */
		  else
			{
			  /* b1. Unknown frame: 3 copies */
			  if (three)
				{
				  /* Replication: 3 SRs for the same coordinates */
				  /* Frames 3, 4 and 5 equal 0, 1 and 2 (reverse) */
				  for (frame=FRAMES; frame < 2*FRAMES; frame++)
					{
					  /* Setting SR data */
					  sr->sRegions[frame][sr->nRegions[frame]].Pos1 =
                        LengthSequence - pos2 + 1;
					  sr->sRegions[frame][sr->nRegions[frame]].Pos2 = 
                        LengthSequence - pos1 + 1;
					  sr->sRegions[frame][sr->nRegions[frame]].Score = score;

					  sr->nRegions[frame]++;
					}
				}
			  /* b2. Known frame */
			  else
				{
				  /* Convert blast frame 3 into frame 0 to store it */
				  /* Frames 3, 4 and 5 equal 0, 1 and 2 (reverse) */
				  frame = (frame % 3);
				  frame = FRAMES + frame;

				  /* Setting SR data */
				  sr->sRegions[frame][sr->nRegions[frame]].Pos1 =
					LengthSequence - pos2 + 1;
				  sr->sRegions[frame][sr->nRegions[frame]].Pos2 = 
					LengthSequence - pos1 + 1;
				  sr->sRegions[frame][sr->nRegions[frame]].Score = score;

				  sr->nRegions[frame]++;
				}
			} /* end RVS */
  
		  /* ready for next SR... */
		  i = (three)? i + 3: i + 1;
		  three = 0;

		  if ((i + FRAMES) > MAXEVIDENCES)
			printError("Too many SRs: increase MAXSR definition");

		} /* end of if-comment */
    } /* end of while*/

  fclose(file);
  
  return(i);
}
