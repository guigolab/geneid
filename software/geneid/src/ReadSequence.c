/*************************************************************************
*                                                                        *
*   Module: ReadSequence                                                 *
*                                                                        *
*   Reading of sequences of DNA in Fasta format                          *
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

/*  $Id: ReadSequence.c,v 1.1 2000-07-05 08:24:30 eblanco Exp $  */

#include "geneid.h"

extern int VRB;

/* It reads the header of the first DNA sequence (Name) */
int IniReadSequence(FILE* seqfile, char* line)
{
  int res;
  char sAux[MAXSTRING];
  char cAux;

  /* Get locus name */
  res = fscanf(seqfile,"%s",sAux);
  /* Delete '>' */
  strcpy(line,sAux+1);

  /* Jumping until \n of the first fasta line */
  res = fscanf(seqfile,"%c",&cAux);
  while(cAux != '\n')
    res = fscanf(seqfile,"%c",&cAux);
    
  if (res == EOF)
    printError("Bad format DNA sequence\n");  
  
  return(res);
}

/* It reads content of current sequence and the header of next */
int ReadSequence (FILE* seqfile, char* Sequence, char* nextLocus)
{
  long pos;
  int res;
  char mess[MAXSTRING];
  char cAux;

  /* fasta format = "atcgata...atta\n" */
  pos = 0;
  res = fscanf(seqfile,"%s\n",Sequence);
  while((res != EOF) && (Sequence[pos] != '>'))
    { 
      pos = pos + strlen(Sequence + pos);
      res = fscanf(seqfile,"%s",Sequence + pos);

      if ( VRB && !(pos % 100000) )
	{
	  sprintf(mess, "...%ld bp",pos);
	  printMess(mess);
	}
    }

  /* Next Sequence */
  if (res != EOF)
    {
      /* Delete '>' */
      strcpy(nextLocus,Sequence+pos+1);
      Sequence[pos] = '\0';

      /* Jumping until \n of the first fasta line */
      res = fscanf(seqfile,"%c",&cAux);
      while(cAux != '\n')
	res = fscanf(seqfile,"%c",&cAux);
      
      if (res == EOF)
	printError("Bad format DNA sequence\n");  
    }
  return(res);
}

