/*************************************************************************
*                                                                        *
*   Module: PrintSites                                                   *
*                                                                        *
*   Print predicted signals (motif, score and position)                  *
*                                                                        *
*   This file is part of the geneid 1.3 distribution                     *
*                                                                        *
*     Copyright (C) 2003 - Enrique BLANCO GARCIA                         *
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

/*  $Id: PrintSites.c,v 1.6 2003-11-05 14:45:55 eblanco Exp $  */

#include "geneid.h"

extern int GFF;

/* Print a individual signal according to the selected format */
void PrintSite(site* s, int type,
               char Name[], int Strand,
               char* seq, profile* p)
{
  int k;
  char Type[MAXTYPE];
  char strand;
  int offset;
  char sAux[MAXSTRING];
  int i;
  
  /* Corrections because of geneid works with coding (useful) positions */
  /* position + k means "real beginning of the signal" */
  /* position + offset means "real end of the signal" */
  Type[0] = '\0';
  switch(type)
    {
    case ACC: strcpy(Type,sACC);
      k = -1;
      offset = 1;
      break;
    case DON:strcpy(Type,sDON);
      k = 0;      
      offset = 1;
      break;
    case STO:strcpy(Type,sSTO);
      k = 1;
      offset = 2;
      break;
    case STA:strcpy(Type,sSTA);
      k = 0;
      offset = 2;
      break;
	default:
	  /* It will be supposed to be an acceptor splice site */
	  strcpy(Type,sACC);
      k = -1;
      offset = 1;
      break;
    }
  
  /* Depending on the strand: exchange k and offset */
  /* Incomplete signals are filled in with "*" */
  switch(Strand)
    {
    case FORWARD: strand = '+';
      for(i=0; i < p->dimension; i++)
		if (s->Position - p->offset + i >= 0)
		  sAux[i] = seq[s->Position - p->offset + i];
		else
		  sAux[i] = '*';
      sAux[i] = '\0';
      break; 
    case REVERSE: strand = '-';
      k = -k;
      offset = -offset;
      ReverseSubSequence(s->Position - (p->dimension - p->offset - 1), 
						 s->Position + (p->offset), seq, sAux);    
      sAux[p->dimension] = '\0'; 
      break;
	default: 
	  /* It will supposed to be forward */
	  strand = '+';
      for(i=0; i < p->dimension; i++)
		if (s->Position - p->offset + i >= 0)
		  sAux[i] = seq[s->Position - p->offset + i];
		else
		  sAux[i] = '*';
      sAux[i] = '\0';
      break;
    }
  
  /* Effective output for this signal */
  if (GFF)
	{
	  /* Print site: gff format */
	  printf("%s\t%s\t%s\t%ld\t%ld\t%5.2f\t%c\t.\t# %s\n",
			 Name,
			 SITES,
			 Type,
			 (offset>0)? 
			 s->Position+k+COFFSET : 
			 s->Position+k+offset+COFFSET,
			 
			 (offset>0)? 
			 s->Position+k+offset+COFFSET : 
			 s->Position+k+COFFSET,
			 
			 s->Score,
			 strand,
			 sAux);
	}
  else
	{
	  /* Print site: default format */
	  printf("%8s %8ld %8ld\t%5.2f\t%c\t%s\n",
			 Type,
			 
			 (offset>0)? 
			 s->Position+k+COFFSET : 
			 s->Position+k+offset+COFFSET,
			 
			 (offset>0)? 
			 s->Position+k+offset+COFFSET : 
			 s->Position+k+COFFSET,
			 
			 s->Score,
			 strand,
			 sAux);
	}
}

/* Print a list of signals according to the selected format */
void PrintSites (site* s,
                 long ns,
                 int type,
                 char Name[],
                 int Strand,
                 long l1, long l2,
                 char* seq,
                 profile *p)
{
  long i;
  char Type[MAXTYPE];
  char strand;
  long limit;
  long printPoint;
  
  /* Skip repeated sites because of overlapping fragments */
  /* depending on the type of the output signals */
  limit = l1 + OVERLAP;
  i = 0;
  printPoint = ns;
  Type[0] = '\0';
  switch(type)
    {
    case ACC: strcpy(Type,sACC);
      if ((Strand == REVERSE) && (l1))
		{
		  i = ns-1;
		  while(i>=0 && ((s+i)->Position < limit))
			i--;
		  printPoint = i+1;
		  i = 0;
		}
      break;
	  
    case STA:strcpy(Type,sSTA);
      if ((Strand == REVERSE) && (l1))
		{
		  i = ns-1;
		  while(i>=0 && ((s+i)->Position < limit))
			i--;
		  printPoint = i+1;
		  i = 0;
		}
      break;
	  
    case DON:strcpy(Type,sDON);
      if ((Strand == FORWARD) && (l1))
		{
		  i = 0;
		  while(i<ns && ((s+i)->Position < limit))
			i++;
		}
      break;
      
    case STO:strcpy(Type,sSTO);
      if ((Strand == FORWARD) && (l1))
		{
		  i = 0;
		  while(i<ns && ((s+i)->Position < limit))
			i++;
		}
      if ((Strand == REVERSE) && (l1))
		{
		  i = ns-1;
		  while(i>=0 && ((s+i)->Position < limit))
			i--;
		  printPoint = i+1;
		  i = 0;
		}
      break;
    }
  
  strand = (Strand == FORWARD)? '+' : '-';
  
  /* Output header for the input (list, fragment) */
  printf("# %ss(%c) predicted in sequence %s: [%ld,%ld]\n",
		 Type,
		 strand,
		 Name, l1, l2);
  
  /* Print part of the list by computing the values i and printPoint */
  for (; i < printPoint; i++)
    PrintSite((s+i),type,Name,Strand,seq,p);
}


