/*************************************************************************
*                                                                        *
*   Module: ReadGeneModel                                                *
*                                                                        *
*   Reading of the rules used at gene construction.                      *
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

/*  $Id: ReadGeneModel.c,v 1.1 2000-07-05 08:22:42 eblanco Exp $  */

#include "geneid.h"

long ReadGeneModel (FILE *file, dict *d, int nc[], int ne[], 
		    int UC[][MAXENTRY], int DE[][MAXENTRY], 
          long md[], long Md[], int block[])
{
  char line[MAXLINE];
  char *t1;
  char *line1;
  char *line2;
  char *line3;
  char *line4;
  int a;
  int nlines=0;

  /* GeneModel format is: 
     F1:F2:(...):Fn   F1:F2:(...):Fm dmin:dMax  [block|...] */
  
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      /* For each line extract the features (upstream/downstream), */
      /* the minMax distances and the (optional) blocks */
      if (line[0] !='#')
	{
	  /* Split line in four parts: UC DE Dist block */
	  line1 = (char *) strtok(line," ");
	  line2 = (char *) strtok(NULL," ");
	  line3 = (char *) strtok(NULL," "); 
	  line4 = (char *) strtok(NULL," ");

	  if (line1 == NULL || line2 == NULL || line3 == NULL)
	    printError("Bad format: GeneModel rule");

	  /* For each feature in this upstream set */
	  for ( t1 =(char *) strtok(line1,":");
		t1 != NULL;
		t1 = (char *) strtok(NULL, ":") )
	    { 
	      a = setkeyDict(d,t1);
	      UC[a][nc[a]++] = nlines;
	    }
	  
	  /* For each feature in this downstream set */
	  for ( t1 =(char *) strtok(line2,":");
		t1 != NULL;
		t1 = (char *) strtok(NULL, ":") )
	    { 
	      a = setkeyDict(d,t1); 
	      DE[a][ne[a]++]=nlines;
	    } 
	  
	  /* Read the distances */
	  t1 =(char *) strtok(line3,":");
	  
	  if (t1 == NULL)
	    printError("Bad format: distances rule");

	  md[nlines] = atol(t1); 
	  t1 = (char *) strtok(NULL, ":");

	  if (t1 == NULL)
	    printError("Bad format: distances rule");

	  Md[nlines] = atol(t1); 
	  
	  /* Read the block ... if exists */
	  if (line4 != NULL) 
	    {
         block[nlines] = BLOCK;
	    }
	  else 
	    block[nlines] = NONBLOCK;
	  nlines++;
	}
    }
  return(nlines);
}


